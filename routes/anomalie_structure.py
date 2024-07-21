from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationAnomalieStructure
from utils import role_requis

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
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Anomalie Structure analysis...\"\n"
        script_content += f"sniffles --input \"{config['input_dir']}\" --vcf \"{config['output_dir']}\"\n"
        script_content += f"echo \"Anomalie Structure analysis completed.\"\n"
    return jsonify(script=script_content)

@anomalie_structure_bp.route('/download_anomalie_structure_script', methods=['GET'])
@role_requis('superadmin')
def download_anomalie_structure_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_anomalie_structure:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Anomalie Structure analysis...\"\n"
        script_content += f"sniffles --input \"{config['input_dir']}\" --vcf \"{config['output_dir']}\"\n"
        script_content += f"echo \"Anomalie Structure analysis completed.\"\n"
        
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

