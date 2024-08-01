import subprocess
import shlex
import os
from flask import Blueprint, render_template, request, jsonify, make_response, send_file
from extensions import db
from models import ConfigurationMerge, Workflow
from utils import role_requis
from datetime import datetime
import json

merge_bp = Blueprint('merge_bp', __name__)

configurations_merge = []

@merge_bp.route('/bam_merger', methods=['GET', 'POST'])
@role_requis('superadmin')
def bam_merger():
    if request.method == 'POST':
        input_dir = request.form['input_dir']
        output_dir = request.form['output_dir']
        
        if not all([input_dir, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_merge.append({
            "input_dir": input_dir,
            "output_dir": output_dir
        })
        
        configurations_merge_db = ConfigurationMerge(
            input_dir=input_dir,
            output_dir=output_dir
        )
        db.session.add(configurations_merge_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('bam_merger.html')

@merge_bp.route('/generate_bam_script', methods=['GET'])
@role_requis('superadmin')
def generate_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        merge_dir = f"{config['output_dir']}/merge-BAM"
        log_file = f"{merge_dir}/merge_log.txt"
        report_file = f"{merge_dir}/merge_report.html"
        status_file = f"{merge_dir}/merge_status.txt"
        merged_bam = f"{merge_dir}/merge.bam"

        script_content += f"mkdir -p \"{merge_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting merge for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"samtools merge \"{merged_bam}\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Merging complete for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Merging failed for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Log Report</title></head><body><div class=\"log-container\"><h1>Log Report</h1>' > {report_file}\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> {report_file}\n"
        script_content += f"done < {log_file}\n"
        script_content += f"echo '</div></body></html>' >> {report_file}\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)



@merge_bp.route('/download_bam_script', methods=['GET'])
@role_requis('superadmin')
def download_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        merge_dir = f"{config['output_dir']}/merge-BAM"
        log_file = f"{merge_dir}/merge_log.txt"
        report_file = f"{merge_dir}/merge_report.html"
        status_file = f"{merge_dir}/merge_status.txt"
        merged_bam = f"{merge_dir}/merge.bam"

        script_content += f"mkdir -p \"{merge_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting merge for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"samtools merge \"{merged_bam}\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Merging complete for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Merging failed for BAM files in {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Log Report</title></head><body><div class=\"log-container\"><h1>Log Report</h1>' > {report_file}\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> {report_file}\n"
        script_content += f"done < {log_file}\n"
        script_content += f"echo '</div></body></html>' >> {report_file}\n"
    
    script_path = '/data/Script_Site/tmp/bam_merge_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="bam_merge_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=bam_merge_script.sh"
    return response



@merge_bp.route('/get_configurations_merge', methods=['GET'])
@role_requis('superadmin')
def get_configurations_merge():
    return jsonify(configurations_merge)

@merge_bp.route('/delete_config_merge', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_merge():
    index = request.json['index']
    try:
        configurations_merge.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration does not exist.")

@merge_bp.route('/start_merge_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="BAM Merge", status="Running", start_time=datetime.utcnow(), output_dir=configurations_merge[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/bam_merge_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_merge[-1]['output_dir'] + "/merge_status.txt"
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

@merge_bp.route('/history-merge')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationMerge.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-merge.html', configurations=configurations)