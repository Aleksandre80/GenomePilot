from flask import Blueprint, render_template, request, jsonify, make_response, send_file
from extensions import db
from models import ConfigurationReadsLength
from utils import role_requis

reads_length_bp = Blueprint('reads_length_bp', __name__)

configurations_reads_length = []

@reads_length_bp.route('/reads_length', methods=['GET','POST'])
@role_requis('superadmin')
def reads_length():
    if request.method == 'POST':
        input_file = request.form.get('input_file')
        output_dir = request.form.get('output_dir')
        length_option = request.form.get('length_option')
        min_length = request.form.get('min_length')
        min_length_between = request.form.get('min_length_between')
        max_length_between = request.form.get('max_length_between')

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
    configurations = ConfigurationReadsLength.query.all()
    configurations_list = [{
        "id": config.id,
        "input_file": config.input_file,
        "output_dir": config.output_dir,
        "length_option": config.length_option,
        "min_length": config.min_length,
        "max_length": config.max_length,
        "date_created": config.date_created.isoformat()
    } for config in configurations]
    
    return jsonify(configurations=configurations_list), 200

@reads_length_bp.route('/delete_config_reads_length', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_reads_length():
    index = request.json['index']
    try:
        configurations_reads_length.pop(index)
        return jsonify(success=True, configurations=configurations_reads_length)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")

# Routes pour générer et télécharger les scripts
@reads_length_bp.route('/generate_reads_length_script', methods=['GET'])
@role_requis('superadmin')
def generate_reads_length_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_reads_length:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"python bam_quality_filter.py \"{config['input_file']}\" \"{config['output_dir']}\" \"{config['chr']}\" \"{config['pos1']}\" \"{config['pos2']}\" \"{config['phred_min']}\" \"{config['logs']}/logs_phred_filter.log\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
    return jsonify(script=script_content)

@reads_length_bp.route('/download_reads_length_script', methods=['GET'])
@role_requis('superadmin')
def download_reads_length_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_reads_length:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"echo \"Starting Calculing reads average...\"\n"
        script_content += f"python bam_quality_filter.py \"{config['input_file']}\" \"{config['output_dir']}\" \"{config['chr']}\" \"{config['pos1']}\" \"{config['pos2']}\" \"{config['phred_min']}\" \"{config['logs']}/logs_phred_filter.log\"\n"
        script_content += f"echo \"Calcul reads average completed.\"\n"
        
    script_path = '/data/Script_Site/tmp/bam_quality_filter.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="bam_quality_filter.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=bam_quality_filter.sh"
    return response
