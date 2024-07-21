from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationMerge
from utils import role_requis

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
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"samtools merge \"{config['output_dir']}/merged.bam\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"echo \"Merging complete for BAM files in {config['input_dir']}\"\n\n"
    return jsonify(script=script_content)

@merge_bp.route('/download_bam_script', methods=['GET'])
@role_requis('superadmin')
def download_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"samtools merge \"{config['output_dir']}/merged.bam\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"echo \"Merging complete for BAM files in {config['input_dir']}\"\n\n"
    
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
