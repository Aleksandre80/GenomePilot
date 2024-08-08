from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationMethylation, Workflow
from utils import role_requis
import json
import os
import subprocess
import shlex
from datetime import datetime
import re

methylation_bp = Blueprint('methylation_bp', __name__)

configurations_methylation = []

@methylation_bp.route('/methylation', methods=['GET', 'POST'])
@role_requis('superadmin')
def methylation():
    if request.method == 'POST':
        input_dir = request.form['input_dir']
        output_dir = request.form['output_dir']
        ref_genome = request.form['ref_genome']
        methylationModelMethyl = request.form['methylationModel']
        
        # Calculate the basic model from the selected methylation model
        methylationModelBasic = re.match(r"^(.*?@v\d+\.\d+\.\d+)", methylationModelMethyl).group(1)
        
        if not all([input_dir, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_methylation.append({
            "input_dir": input_dir,
            "output_dir": output_dir,
            "ref_genome" : ref_genome,
            "methylationModelBasic" : methylationModelBasic,
            "methylationModelMethyl" : methylationModelMethyl
        })
        
        configurations_methylation_db = ConfigurationMethylation(
            input_dir=input_dir,
            output_dir=output_dir,
            ref_genome = ref_genome,
            methylationModelBasic = methylationModelBasic,
            methylationModelMethyl = methylationModelMethyl
        )
        db.session.add(configurations_methylation_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('methylation.html')

@methylation_bp.route('/generate_methylation_script', methods=['GET'])
@role_requis('superadmin')
def generate_methylation_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    
    for config in configurations_methylation:
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        model_basename = os.path.basename(config['methylationModelBasic']).replace('.bin', '')
        
        output_dir = os.path.join(config['output_dir'], "Methylation")
        log_file = f"{output_dir}/methylation_log.txt"
        report_file = f"{output_dir}/methylation_report.html"
        status_file = f"{output_dir}/methylation_status.txt"
        output_bam = f"{output_dir}/{ref_basename}_{model_basename}.bam"
        
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting methylation analysis for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter la méthylation
        script_content += f"/home/grid/dorado-0.7.2-linux-x64/bin/dorado basecaller \\\n"
        script_content += f"    \"/home/grid/dorado-0.7.2-linux-x64/bin/{config['methylationModelBasic']}\" \\\n"
        script_content += f"    \"{config['input_dir']}\" \\\n"
        script_content += f"    --modified-bases-models \"/home/grid/dorado-0.7.2-linux-x64/bin/{config['methylationModelMethyl']}\" \\\n"
        script_content += f"    --reference \"{config['ref_genome']}\" \\\n"
        script_content += f"    > \"{output_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Methylation analysis completed successfully.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Methylation analysis failed.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Methylation Log Report</title></head><body><div class=\"log-container\"><h1>Methylation Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)
    
    return jsonify(script=escaped_script_content)





@methylation_bp.route('/download_methylation_script', methods=['GET'])
@role_requis('superadmin')
def download_methylation_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_methylation:
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        model_basename = os.path.basename(config['methylationModelBasic']).replace('.bin', '')
        
        output_dir = os.path.join(config['output_dir'], "Methylation")
        log_file = f"{output_dir}/methylation_log.txt"
        report_file = f"{output_dir}/methylation_report.html"
        status_file = f"{output_dir}/methylation_status.txt"
        output_bam = f"{output_dir}/{ref_basename}_{model_basename}.bam"
        
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting methylation analysis for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter la méthylation
        script_content += f"/home/grid/dorado-0.7.2-linux-x64/bin/dorado basecaller \\\n"
        script_content += f"    \"/home/grid/dorado-0.7.2-linux-x64/bin/{config['methylationModelBasic']}\" \\\n"
        script_content += f"    \"{config['input_dir']}\" \\\n"
        script_content += f"    --modified-bases-models \"/home/grid/dorado-0.7.2-linux-x64/bin/{config['methylationModelMethyl']}\" \\\n"
        script_content += f"    --reference \"{config['ref_genome']}\" \\\n"
        script_content += f"    > \"{output_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Methylation analysis completed successfully.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Methylation analysis failed.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Methylation Log Report</title></head><body><div class=\"log-container\"><h1>Methylation Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
        
    script_path = '/data/Script_Site/tmp/methylation_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="methylation_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=methylation_script.sh"
    return response

@methylation_bp.route('/get_configurations_methylation', methods=['GET'])
@role_requis('superadmin')
def get_configurations_methylation():
    return jsonify(configurations_methylation)

@methylation_bp.route('/delete_config_methylation', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_methylation():
    index = request.json['index']
    try:
        configurations_methylation.pop(index)
        return jsonify(success=True, configurations=configurations_methylation)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

@methylation_bp.route('/start_methylation_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_methylation_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Basecalling Methylation", status="Running", start_time=datetime.utcnow(), output_dir=configurations_methylation[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/methylation_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_methylation[-1]['output_dir'] + "/Methylation/methylation_status.txt"
            if os.path.exists(status_file):
                with open(status_file, 'r') as file:
                    status_info = file.read().strip()
                    status, end_time = status_info.split(' - ')
                    new_workflow.status = "Completed" if status == "completed" else "Failed"
                    new_workflow.end_time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
            else:
                new_workflow.status = "Failed"
                new_workflow.end_time = datetime.utcnow()
            
            db.session.commit()

        except Exception as e:
            print(f"Error: {e}")
            workflow = Workflow.query.get(new_workflow.id)
            workflow.status = "Failed"
            workflow.end_time = datetime.utcnow()
            db.session.commit()

        return jsonify(success=True, report=new_workflow.status)
    
    return jsonify(success=False, message="Invalid request method. Use POST.")
