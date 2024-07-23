from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationAnomalieStructure, Workflow
from utils import role_requis
import json
import subprocess
import shlex

anomalie_structure_bp = Blueprint('anomalie_structure_bp', __name__)

configurations_anomalie_structure = []

@anomalie_structure_bp.route('/anomalie_structure', methods=['GET', 'POST'])
@role_requis('superadmin')
def anomalie_structure():
    if request.method == 'POST':
        input_dir = request.form['input_dir']
        output_dir = request.form['output_dir']
        
        if not all([input_dir, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_anomalie_structure.append({
            "input_dir": input_dir,
            "output_dir": output_dir
        })
        
        configurations_anomalie_structure_db = ConfigurationAnomalieStructure(
            input_dir=input_dir,
            output_dir=output_dir
        )
        db.session.add(configurations_anomalie_structure_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('anomalie-structure.html')

@anomalie_structure_bp.route('/generate_anomalie_structure_script', methods=['GET'])
@role_requis('superadmin')
def generate_anomalie_structure_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_anomalie_structure:
        log_file = f"{config['output_dir']}/anomalie_structure_log.txt"
        report_file = f"{config['output_dir']}/anomalie_structure_report.html"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting Anomalie Structure analysis for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"source ~/.pyenv/versions/sniffles-env/bin/activate\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"cd /usr/local/bin/Sniffles-2.4\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Environment activated and directory changed to /usr/local/bin/Sniffles-2.4\" >> \"{log_file}\"\n"
        
        # Commande pour générer le fichier VCF
        script_content += f"./sniffles --input \"{config['input_dir']}\" --vcf \"{config['output_dir']}/output.vcf\" >> \"{log_file}\" 2>&1\n"
        
        script_content += f"if [ -f \"{config['output_dir']}/output.vcf\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Anomalie Structure analysis completed and output.vcf generated.\" >> \"{log_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Anomalie Structure analysis failed. output.vcf not generated.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Anomalie Structure Log Report</title></head><body><div class=\"log-container\"><h1>Anomalie Structure Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)




@anomalie_structure_bp.route('/download_anomalie_structure_script', methods=['GET'])
@role_requis('superadmin')
def download_anomalie_structure_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_anomalie_structure:
        log_file = f"{config['output_dir']}/anomalie_structure_log.txt"
        report_file = f"{config['output_dir']}/anomalie_structure_report.html"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting Anomalie Structure analysis for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        script_content += f"source ~/.pyenv/versions/sniffles-env/bin/activate\n"
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"cd /usr/local/bin/Sniffles-2.4\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Environment activated and directory changed to /usr/local/bin/Sniffles-2.4\" >> \"{log_file}\"\n"
        
        # Commande pour générer le fichier VCF
        script_content += f"./sniffles --input \"{config['input_dir']}\" --vcf \"{config['output_dir']}/output.vcf\" >> \"{log_file}\" 2>&1\n"
        
        script_content += f"if [ -f \"{config['output_dir']}/output.vcf\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Anomalie Structure analysis completed and output.vcf generated.\" >> \"{log_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Anomalie Structure analysis failed. output.vcf not generated.\" >> \"{log_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Anomalie Structure Log Report</title></head><body><div class=\"log-container\"><h1>Anomalie Structure Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
        
    script_path = '/data/Script_Site/tmp/anomalie_structure_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="anomalie_structure_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=anomalie_structure_script.sh"
    return response

@anomalie_structure_bp.route('/get_configurations_anomalie_structure', methods=['GET'])
@role_requis('superadmin')
def get_configurations_anomalie_structure():
    return jsonify(configurations_anomalie_structure)

@anomalie_structure_bp.route('/delete_config_anomalie_structure', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_anomalie_structure():
    index = request.json['index']
    try:
        configurations_anomalie_structure.pop(index)
        return jsonify(success=True, configurations=configurations_anomalie_structure)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

@anomalie_structure_bp.route('/start_anomalie_structure_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Anomalie de Structure", status="Running")
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/anomalie_structure_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            # Assuming report_file path is stored in a way that it can be dynamically resolved
            report_file = configurations_anomalie_structure[-1]['output_dir'] + "/anomalie_structure_report.html"
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