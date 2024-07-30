from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationMethylartistViolon, Workflow
from utils import role_requis
import json
import os
import subprocess
import shlex
from datetime import datetime

methylartist_violon_bp = Blueprint('methylartist_violon_bp', __name__)

configurations_methylartist_violon = []

@methylartist_violon_bp.route('/methylartist_violon', methods=['GET', 'POST'])
@role_requis('superadmin')
def methylartist_violon():
    if request.method == 'POST':
        input_file = request.form['input_file']
        output_dir = request.form['output_dir']
        ref_genome = request.form['ref_genome']
        ref_bed = request.form['ref_bed']
        
        
        if not all([input_file, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_methylartist_violon.append({
            "input_file": input_file,
            "output_dir": output_dir,
            "ref_genome": ref_genome,
            "ref_bed" : ref_bed
        })
        
        configurations_methylartist_violon_db = ConfigurationMethylartistViolon(
            input_file=input_file,
            output_dir=output_dir,
            ref_genome=ref_genome,
            ref_bed = ref_bed
        )
        db.session.add(configurations_methylartist_violon_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('methylartist_violon.html')

@methylartist_violon_bp.route('/generate_methylartist_violon_script', methods=['GET'])
@role_requis('superadmin')
def generate_methylartist_violon_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate methplotlib_env\n\n"
    
    for config in configurations_methylartist_violon:
        log_file = f"{config['output_dir']}/methylartist_violon_log.txt"
        report_file = f"{config['output_dir']}/methylartist_violon_report.html"
        input_bam = config['input_file']  # Using input_file directly from the config
        ref_bed = config['ref_bed']
        ref_genome = config['ref_genome']
        
        bam_basename = os.path.basename(input_bam).replace('.bam', '')
        bed_basename = os.path.basename(ref_bed).replace('.bed', '')
        output_tsv = f"{config['output_dir']}/{bed_basename}.{bam_basename}.tsv"
        output_png = f"{config['output_dir']}/{bed_basename}.{bam_basename}-violon.svg"
        
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting methylartist_violon analysis for input file {input_bam}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter methylartist_violon segmeth
        script_content += f"methylartist_violon segmeth -b \"{input_bam}\" -i \"{ref_bed}\" -p 32 --ref \"{ref_genome}\" --motif CG > \"{output_tsv}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - methylartist_violon segmeth completed successfully.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter segplot
        script_content += f"    segplot -s \"{output_tsv}\" -v -o \"{output_png}\" 2>> \"{log_file}\"\n"
        
        script_content += f"    if [ $? -eq 0 ]; then\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - segplot completed successfully.\" >> \"{log_file}\"\n"
        script_content += f"    else\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - segplot failed.\" >> \"{log_file}\"\n"
        script_content += f"    fi\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - methylartist_violon segmeth failed.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>methylartist_violon Log Report</title></head><body><div class=\"log-container\"><h1>methylartist_violon Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)
    
    return jsonify(script=escaped_script_content)






@methylartist_violon_bp.route('/download_methylartist_violon_script', methods=['GET'])
@role_requis('superadmin')
def download_methylartist_violon_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate methplotlib_env\n\n"
    for config in configurations_methylartist_violon:
        log_file = f"{config['output_dir']}/methylartist_violon_log.txt"
        report_file = f"{config['output_dir']}/methylartist_violon_report.html"
        input_bam = config['input_file']  # Using input_file directly from the config
        ref_bed = config['ref_bed']
        ref_genome = config['ref_genome']
        
        bam_basename = os.path.basename(input_bam).replace('.bam', '')
        bed_basename = os.path.basename(ref_bed).replace('.bed', '')
        output_tsv = f"{config['output_dir']}/{bed_basename}.{bam_basename}.tsv"
        output_png = f"{config['output_dir']}/{bed_basename}.{bam_basename}-violon.svg"
        
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting methylartist_violon analysis for input file {input_bam}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter methylartist_violon segmeth
        script_content += f"methylartist_violon segmeth -b \"{input_bam}\" -i \"{ref_bed}\" -p 32 --ref \"{ref_genome}\" --motif CG > \"{output_tsv}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - methylartist_violon segmeth completed successfully.\" >> \"{log_file}\"\n"
        
        # Commande pour exécuter segplot
        script_content += f"    segplot -s \"{output_tsv}\" -v -o \"{output_png}\" 2>> \"{log_file}\"\n"
        
        script_content += f"    if [ $? -eq 0 ]; then\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - segplot completed successfully.\" >> \"{log_file}\"\n"
        script_content += f"    else\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - segplot failed.\" >> \"{log_file}\"\n"
        script_content += f"    fi\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - methylartist_violon segmeth failed.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>methylartist_violon Log Report</title></head><body><div class=\"log-container\"><h1>methylartist_violon Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    script_path = '/data/Script_Site/tmp/methylartist_violon_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="methylartist_violon_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=methylartist_violon_script.sh"
    return response

@methylartist_violon_bp.route('/get_configurations_methylartist_violon', methods=['GET'])
@role_requis('superadmin')
def get_configurations_methylartist_violon():
    return jsonify(configurations_methylartist_violon)

@methylartist_violon_bp.route('/delete_config_methylartist_violon', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_methylartist_violon():
    index = request.json['index']
    try:
        configurations_methylartist_violon.pop(index)
        return jsonify(success=True, configurations=configurations_methylartist_violon)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

@methylartist_violon_bp.route('/start_methylartist_violon_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_methylartist_violon_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="methylartist_violon", status="Running", start_time=datetime.utcnow(), output_dir=configurations_methylartist_violon[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/methylartist_violon_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_methylartist_violon[-1]['output_dir'] + "/methylartist_status.txt"
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
