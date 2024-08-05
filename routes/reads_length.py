from flask import Blueprint, render_template, request, jsonify, make_response, send_file
from extensions import db
from models import ConfigurationReadsLength, Workflow
from utils import role_requis
import json
import subprocess
import shlex
import os
from datetime import datetime

reads_length_bp = Blueprint('reads_length_bp', __name__)

configurations_reads_length = []

@reads_length_bp.route('/reads_length', methods=['GET', 'POST'])
@role_requis('superadmin')
def reads_length():
    if request.method == 'POST':
        input_file = request.form.get('input_file')
        output_dir = request.form.get('output_dir')
        length_option = request.form.get('length_option')
        min_length = request.form.get('min_length')
        min_length_between = request.form.get('min_length_between')
        max_length_between = request.form.get('max_length_between')

        # Vérifiez que les champs nécessaires sont définis en fonction de l'option sélectionnée
        if length_option == 'sup' and not min_length:
            return jsonify(success=False, message="Veuillez fournir une longueur minimale pour l'option 'sup'."), 400
        elif length_option == 'between' and (not min_length_between or not max_length_between):
            return jsonify(success=False, message="Veuillez fournir des longueurs minimale et maximale pour l'option 'between'."), 400

        new_config = ConfigurationReadsLength(
            input_file=input_file,
            output_dir=output_dir,
            length_option=length_option,
            min_length=min_length if length_option == 'sup' else min_length_between,
            max_length=max_length_between if length_option == 'between' else None
        )
        
        db.session.add(new_config)
        db.session.commit()
        
        configurations_reads_length.append({
            "id": new_config.id,
            "input_file": new_config.input_file,
            "output_dir": new_config.output_dir,
            "length_option": new_config.length_option,
            "min_length": new_config.min_length,
            "max_length": new_config.max_length,
            "date_created": new_config.date_created.isoformat()
        })
        
        return jsonify(success=True, configurations=configurations_reads_length), 200
    return render_template('reads-length.html')


@reads_length_bp.route('/get_configurations_reads_length', methods=['GET'])
@role_requis('superadmin')
def get_configurations_reads_length():
    return jsonify(configurations_reads_length)
    
@reads_length_bp.route('/delete_config_reads_length', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_reads_length():
    index = request.json['index']
    try:
        configurations_reads_length.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration not found.")

@reads_length_bp.route('/generate_reads_length_script', methods=['GET'])
@role_requis('superadmin')
def generate_reads_length_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_reads_length:
        bam_filename = os.path.basename(config['input_file'])
        bam_basename = bam_filename.replace('.bam', '')
        
        if config['length_option'] == 'sup':
            length_desc = f"sup_{config['min_length']}"
        elif config['length_option'] == 'between':
            length_desc = f"between_{config['min_length_between']}_{config['max_length_between']}"
        else:
            length_desc = "unknown_length"

        output_dir = os.path.join(config['output_dir'], "Longueur_Reads")
        log_file = f"{output_dir}/reads_length_log.txt"
        report_file = f"{output_dir}/reads_length_report.html"
        status_file = f"{output_dir}/reads_length_status.txt"
        filtered_bam = f"{output_dir}/{bam_basename}_{length_desc}.bam"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BAM filtering by read length for input file {config['input_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour filtrer le fichier BAM selon l'option de longueur
        if config['length_option'] == 'sup':
            script_content += f"samtools view -h \"{config['input_file']}\" | awk '{{if($1 ~ /^@/ || length($10) >= {config['min_length']}) print}}' | samtools view -b -o \"{filtered_bam}\" 2>> \"{log_file}\"\n"
        elif config['length_option'] == 'between':
            script_content += f"samtools view -h \"{config['input_file']}\" | awk '{{if($1 ~ /^@/ || (length($10) >= {config['min_length_between']} && length($10) <= {config['max_length_between']})) print}}' | samtools view -b -o \"{filtered_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{filtered_bam}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering by read length completed and {bam_basename}_{length_desc}.bam generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering by read length failed. {bam_basename}_{length_desc}.bam not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Reads Length Filtering Log Report</title></head><body><div class=\"log-container\"><h1>Reads Length Filtering Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)





@reads_length_bp.route('/download_reads_length_script', methods=['GET'])
@role_requis('superadmin')
def download_reads_length_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_reads_length:
        bam_filename = os.path.basename(config['input_file'])
        bam_basename = bam_filename.replace('.bam', '')
        
        if config['length_option'] == 'sup':
            length_desc = f"sup_{config['min_length']}"
        elif config['length_option'] == 'between':
            length_desc = f"between_{config['min_length_between']}_{config['max_length_between']}"
        else:
            length_desc = "unknown_length"

        output_dir = os.path.join(config['output_dir'], "Longueur_Reads")
        log_file = f"{output_dir}/reads_length_log.txt"
        report_file = f"{output_dir}/reads_length_report.html"
        status_file = f"{output_dir}/reads_length_status.txt"
        filtered_bam = f"{output_dir}/{bam_basename}_{length_desc}.bam"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting BAM filtering by read length for input file {config['input_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour filtrer le fichier BAM selon l'option de longueur
        if config['length_option'] == 'sup':
            script_content += f"samtools view -h \"{config['input_file']}\" | awk '{{if($1 ~ /^@/ || length($10) >= {config['min_length']}) print}}' | samtools view -b -o \"{filtered_bam}\" 2>> \"{log_file}\"\n"
        elif config['length_option'] == 'between':
            script_content += f"samtools view -h \"{config['input_file']}\" | awk '{{if($1 ~ /^@/ || (length($10) >= {config['min_length_between']} && length($10) <= {config['max_length_between']})) print}}' | samtools view -b -o \"{filtered_bam}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{filtered_bam}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering by read length completed and {bam_basename}_{length_desc}.bam generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - BAM filtering by read length failed. {bam_basename}_{length_desc}.bam not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Reads Length Filtering Log Report</title></head><body><div class=\"log-container\"><h1>Reads Length Filtering Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
        
    script_path = '/data/Script_Site/tmp/reads_length_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="reads_length_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=reads_length_script.sh"
    return response


@reads_length_bp.route('/start_reads_length_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_reads_length_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Reads Length Filtering", status="Running", start_time=datetime.utcnow(), output_dir=configurations_reads_length[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/reads_length_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_reads_length[-1]['output_dir'] + "/Longueur_Reads/reads_length_status.txt"
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

@reads_length_bp.route('/history-reads-length')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationReadsLength.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-reads-length.html', configurations=configurations)