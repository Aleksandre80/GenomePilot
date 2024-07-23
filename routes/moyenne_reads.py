from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationMoyenneReads
from utils import role_requis

moyenne_reads_bp = Blueprint('moyenne_reads_bp', __name__)

configurations_moyenne_reads = []

@moyenne_reads_bp.route('/moyenne_reads', methods=['GET', 'POST'])
@role_requis('superadmin')
def moyenne_reads():
    if request.method == 'POST':
        input_file = request.form['input_file']
        output_dir = request.form['output_dir']
        bed_file = request.form['bed_file']
        
        if not all([input_file, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_moyenne_reads.append({
            "input_dir": input_file,
            "output_dir": output_dir,
            "bed_file": bed_file
        })
        
        configurations_moyenne_reads_db = ConfigurationMoyenneReads(
            input_file=input_file,
            output_dir=output_dir,
            bed_file=bed_file
        )
        db.session.add(configurations_moyenne_reads_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('moyenne-reads.html')

@moyenne_reads_bp.route('/generate_moyenne_reads_script', methods=['GET'])
@role_requis('superadmin')
def generate_anomalie_structure_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_moyenne_reads:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"bedtools coverage -a \"{config['bed_file']}\" -b \"{config['input_dir']}\" > \"{config['output_dir']}/coverage.txt\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
    return jsonify(script=script_content)

@moyenne_reads_bp.route('/download_moyenne_reads_script', methods=['GET'])
@role_requis('superadmin')
def download_moyenne_reads_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_moyenne_reads:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"bedtools coverage -a \"{config['bed_file']}\" -b \"{config['input_dir']}\" > \"{config['output_dir']}/coverage.txt\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
        
    script_path = '/data/Script_Site/tmp/moyenne_reads_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="moyenne_reads_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=moyenne_reads_script.sh"
    return response

@moyenne_reads_bp.route('/get_configurations_moyenne_reads', methods=['GET'])
@role_requis('superadmin')
def get_configurations_moyenne_reads():
    return jsonify(configurations_moyenne_reads)

@moyenne_reads_bp.route('/delete_config_moyenne_reads', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_moyenne_reads():
    index = request.json['index']
    try:
        configurations_moyenne_reads.pop(index)
        return jsonify(success=True, configurations=configurations_moyenne_reads)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

