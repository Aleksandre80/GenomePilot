from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationCibleReads, Workflow
from utils import role_requis
import json
import os
import subprocess
import shlex

cible_reads_bp = Blueprint('cible_reads_bp', __name__)

configurations_cible_reads = []

@cible_reads_bp.route('/cible_reads', methods=['GET', 'POST'])
@role_requis('superadmin')
def cible_reads():
    if request.method == 'POST':
        input_file = request.form['input_file']
        output_dir = request.form['output_dir']
        bed_file = request.form['bed_file']
        
        if not all([input_file, bed_file,output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_cible_reads.append({
            "input_file": input_file,
            "output_dir": output_dir,
            "bed_file": bed_file
        })
        print(configurations_cible_reads)
        
        
        configurations_cible_reads_db = ConfigurationCibleReads(
            input_file=input_file,
            output_dir=output_dir,
            bed_file=bed_file
        )
        db.session.add(configurations_cible_reads_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('cible-reads.html')

@cible_reads_bp.route('/generate_cible_reads_script', methods=['GET'])
@role_requis('superadmin')
def generate_cible_reads_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_cible_reads:
        log_file = f"{config['output_dir']}/cible_reads_log.txt"
        report_file = f"{config['output_dir']}/cible_reads_report.html"
        targeted_bam = f"{config['output_dir']}/targeted.bam"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BAM filtering for regions of interest for input file {config['input_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour filtrer le fichier BAM
        script_content += f"samtools view -b -L \"{config['bed_file']}\" \"{config['input_file']}\" > \"{targeted_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{targeted_bam}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering completed and targeted.bam generated.\" >> \"{log_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering failed. targeted.bam not generated.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>BAM Filtering Log Report</title></head><body><div class=\"log-container\"><h1>BAM Filtering Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)


@cible_reads_bp.route('/download_cible_reads_script', methods=['GET'])
@role_requis('superadmin')
def download_cible_reads_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_cible_reads:
        log_file = f"{config['output_dir']}/cible_reads_log.txt"
        report_file = f"{config['output_dir']}/cible_reads_report.html"
        targeted_bam = f"{config['output_dir']}/targeted.bam"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BAM filtering for regions of interest for input file {config['input_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour filtrer le fichier BAM
        script_content += f"samtools view -b -L \"{config['bed_file']}\" \"{config['input_file']}\" > \"{targeted_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{targeted_bam}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering completed and targeted.bam generated.\" >> \"{log_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering failed. targeted.bam not generated.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>BAM Filtering Log Report</title></head><body><div class=\"log-container\"><h1>BAM Filtering Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
        
    script_path = '/data/Script_Site/tmp/cible_reads_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="cible_reads_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=cible_reads_script.sh"
    return response

@cible_reads_bp.route('/get_configurations_cible_reads', methods=['GET'])
@role_requis('superadmin')
def get_configurations_cible_reads():
    return jsonify(configurations_cible_reads)

@cible_reads_bp.route('/delete_config_cible_reads', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_cible_reads():
    index = request.json['index']
    try:
        configurations_cible_reads.pop(index)
        return jsonify(success=True, configurations=configurations_cible_reads)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

@cible_reads_bp.route('/start_cible_reads_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Filtrage BAM sur cible", status="Running")
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/cible_reads_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            # Assuming report_file path is stored in a way that it can be dynamically resolved
            report_file = configurations_cible_reads[-1]['output_dir'] + "/cible_reads_report.html"
            if os.path.exists(report_file):
                new_workflow.status = "Completed"
            else:
                new_workflow.status = "Completed"
            
            db.session.commit()

        except Exception as e:
            print(f"Error: {e}")
            new_workflow.status = "Completed"
            db.session.commit()

        return jsonify(success=True, report=report_file)
    
    return jsonify(success=False, message="Invalid request method. Use POST.")