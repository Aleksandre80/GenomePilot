from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationPhredQuality
from utils import role_requis
import os

phred_quality_bp = Blueprint('phred_quality_bp', __name__)

configurations_phred_quality = []

@phred_quality_bp.route('/phred_quality', methods=['GET', 'POST'])
@role_requis('superadmin')
def phred_quality():
    if request.method == 'POST':
        input_file = request.form['input_file']
        output_dir = request.form['output_dir']
        chr = request.form['chr']
        pos1 = request.form['pos1']
        pos2 = request.form['pos2']
        phred_min = request.form['phred_min']
        logs = request.form['logs']
        
        if not all([output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_phred_quality.append({
            "input_file": input_file,
            "output_dir": output_dir,
            "chr": chr,
            "pos1": pos1,
            "pos2": pos2,
            "phred_min": phred_min,
            "logs": logs
        })
        
        configurations_phred_quality_db = ConfigurationPhredQuality(
            input_file=input_file,
            output_dir=output_dir,
            chr=chr,
            pos1=pos1,
            pos2=pos2,
            phred_min=phred_min,
            logs=logs
        )
        db.session.add(configurations_phred_quality_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('phred-quality.html')

@phred_quality_bp.route('/generate_phred_quality_script', methods=['GET'])
@role_requis('superadmin')
def generate_anomalie_structure_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_phred_quality:
        output_dir = os.path.join(config['output_dir'], "Phred_Quality")
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"python bam_quality_filter.py \"{config['input_file']}\" \"{output_dir}\" \"{config['chr']}\" \"{config['pos1']}\" \"{config['pos2']}\" \"{config['phred_min']}\" \"{output_dir}/logs_phred_filter.log\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
        #python bam_quality_filter.py exemple.bam filtre.bam chr1 100000 100500 30 log_qualite.log
    return jsonify(script=script_content)

@phred_quality_bp.route('/download_phred_quality_script', methods=['GET'])
@role_requis('superadmin')
def download_phred_quality_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_phred_quality:
        output_dir = os.path.join(config['output_dir'], "Phred_Quality")
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"python bam_quality_filter.py \"{config['input_file']}\" \"{output_dir}\" \"{config['chr']}\" \"{config['pos1']}\" \"{config['pos2']}\" \"{config['phred_min']}\" \"{output_dir}/logs_phred_filter.log\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
        
    script_path = '/data/Script_Site/tmp/bam_quality_filter.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="bam_quality_filter.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=bam_quality_filter.sh"
    return response

@phred_quality_bp.route('/get_configurations_phred_quality', methods=['GET'])
@role_requis('superadmin')
def get_configurations_phred_quality():
    return jsonify(configurations_phred_quality)

@phred_quality_bp.route('/delete_config_phred_quality', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_phred_quality():
    index = request.json['index']
    try:
        configurations_phred_quality.pop(index)
        return jsonify(success=True, configurations=configurations_phred_quality)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

