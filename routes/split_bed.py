from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationSplit, Workflow
import os
from utils import role_requis
import subprocess
import shlex
import json
from datetime import datetime
import re

split_bed_bp = Blueprint('split_bed_bp', __name__)

configurations_split_bed = []

@split_bed_bp.route('/split_bed', methods=['GET', 'POST'])
@role_requis('superadmin') 
def split_bed_creator():
    if request.method == 'POST':
        input_bed = request.form['input_bed']
        output_dir = request.form['output_dir']
        split_size = request.form['split_size']

        configurations_split_bed.append({
            "input_bed": input_bed,
            "output_dir": output_dir,
            "split_size": split_size
        })
        
        configurations_split_bed_db = ConfigurationSplit(
            input_bed=input_bed,
            output_dir=output_dir,
            split_size=split_size
        )
        db.session.add(configurations_split_bed_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('split_bed.html')

@split_bed_bp.route('/generate_split_bed_script', methods=['GET'])
@role_requis('superadmin')
def generate_split_bed_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_split_bed:
        output_dir = os.path.join(config['output_dir'], "Split_Bed")
        input_bed = config['input_bed']
        split_size = config.get('split_size', 1)  # Default to 1 if not provided
        input_bed_basename = os.path.basename(input_bed).replace('.bed', '')
        output_bed = f"{output_dir}/{input_bed_basename}_split_{split_size}.bed"
        log_file = f"{output_dir}/bed_split_log.txt"
        status_file = f"{output_dir}/bed_split_status.txt"
        report_file = f"{output_dir}/bed_split_report.html"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BED splitting for input file {input_bed} with split size {split_size}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour diviser les régions du fichier BED
        script_content += f"awk 'BEGIN {{ OFS = \"\\t\" }} {{ for (i = $2; i < $3; i+={split_size}) {{ print $1, i, (i + {split_size} > $3 ? $3 : i + {split_size}), $4 }} }}' \"{input_bed}\" > \"{output_bed}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BED splitting completed successfully. Output saved to {output_bed}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BED splitting failed.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>BED Split Log Report</title></head><body><div class=\"log-container\"><h1>BED Split Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)


@split_bed_bp.route('/download_split_bed_script', methods=['GET'])
@role_requis('superadmin')
def download_split_bed_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_split_bed:
        output_dir = os.path.join(config['output_dir'], "Split_Bed")
        input_bed = config['input_bed']
        split_size = config.get('split_size', 1)  # Default to 1 if not provided
        input_bed_basename = os.path.basename(input_bed).replace('.bed', '')
        output_bed = f"{output_dir}/{input_bed_basename}_split_{split_size}.bed"
        log_file = f"{output_dir}/bed_split_log.txt"
        status_file = f"{output_dir}/bed_split_status.txt"
        report_file = f"{output_dir}/bed_split_report.html"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BED splitting for input file {input_bed} with split size {split_size}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \{output_dir} \n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour diviser les régions du fichier BED
        script_content += f"awk 'BEGIN {{ OFS = \"\\t\" }} {{ for (i = $2; i < $3; i+={split_size}) {{ print $1, i, (i + {split_size} > $3 ? $3 : i + {split_size}), $4 }} }}' \"{input_bed}\" > \"{output_bed}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BED splitting completed successfully. Output saved to {output_bed}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BED splitting failed.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>BED Split Log Report</title></head><body><div class=\"log-container\"><h1>BED Split Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"


    script_path = '/data/Script_Site/tmp/split_bed_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="split_bed_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=split_bed_script.sh"
    return response


@split_bed_bp.route('/get_configurations_split_bed', methods=['GET'])
@role_requis('superadmin')
def get_configurations_split_bed():
    return jsonify(configurations_split_bed)
    
@split_bed_bp.route('/delete_config_split_bed', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_split_bed():
    index = request.json['index']
    try:
        configurations_split_bed.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration not found.")

@split_bed_bp.route('/start_split_bed_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_split_bed_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Split Bed", status="Running", start_time=datetime.utcnow(), output_dir=configurations_split_bed[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/split_bed_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_split_bed[-1]['output_dir'] + "/Split_Bed/bed_split_status.txt"
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


@split_bed_bp.route('/history-split_bed')
@role_requis('superadmin')
def history():
    configurations = ConfigurationSplit.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-split_bed.html', configurations=configurations)