from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationVCF
import os
from utils import role_requis

vcf_bp = Blueprint('vcf_bp', __name__)

configurations_vcf = []

@vcf_bp.route('/vcf_creator', methods=['GET', 'POST'])
@role_requis('superadmin') 
def vcf_creator():
    if request.method == 'POST':
        ref_genome_path = request.form['ref_genome']
        bam_file = request.form['bam_file']
        output_vcf = request.form['output_vcf']

        configurations_vcf.append({
            "ref_genome": ref_genome_path,
            "bam_file": bam_file,
            "output_vcf": output_vcf
        })
        
        configurations_vcf_db = ConfigurationVCF(
            ref_genome=ref_genome_path,
            bam_file=bam_file,
            output_vcf=output_vcf
        )
        db.session.add(configurations_vcf_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('vcf_creator.html')

@vcf_bp.route('/generate_vcf_script', methods=['GET'])
@role_requis('superadmin') 
def generate_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        output_directory = os.path.dirname(config['output_vcf'])
        vcf_directory = f"{output_directory}/vcf"
        script_content += f"mkdir -p {vcf_directory}\n"
        
        output_vcf_path = f"{vcf_directory}/{os.path.basename(config['output_vcf'])}"
        
        script_content += f"samtools faidx {config['ref_genome']}\n"
        script_content += f"samtools index {config['bam_file']}\n"
        script_content += f"bcftools mpileup -Ou -f {config['ref_genome']} {config['bam_file']} | bcftools call -mv -Ob -o {output_vcf_path}.bcf\n"
        script_content += f"bcftools index {output_vcf_path}.bcf\n"
        script_content += f"bcftools view -Oz -o {output_vcf_path}.vcf.gz {output_vcf_path}.bcf\n"
        script_content += f"tabix -p vcf {output_vcf_path}.vcf.gz\n"
        script_content += f"gunzip -c {output_vcf_path}.vcf.gz > {output_vcf_path}.vcf\n"
        script_content += f"rm -f {output_vcf_path}.bcf {output_vcf_path}.vcf.gz {output_vcf_path}.bcf.csi {output_vcf_path}.vcf.gz.tbi\n"
        script_content += "echo \"Variant calling and file processing completed.\"\n"
    return jsonify(script=script_content)

@vcf_bp.route('/download_vcf_script', methods=['GET'])
@role_requis('superadmin')
def download_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        bam_directory = os.path.dirname(config['bam_file'])
        vcf_directory = os.path.join(bam_directory, "vcf")
        script_content += f"mkdir -p {vcf_directory}\n"

        base_vcf_filename = os.path.basename(config['output_vcf'])
        output_vcf_path = os.path.join(vcf_directory, base_vcf_filename)

        script_content += f"samtools faidx {config['ref_genome']}\n"
        script_content += f"samtools index {config['bam_file']}\n"
        script_content += f"bcftools mpileup -Ou -f {config['ref_genome']} {config['bam_file']} | bcftools call -mv -Ob -o {output_vcf_path}.bcf\n"
        script_content += f"bcftools index {output_vcf_path}.bcf\n"
        script_content += f"bcftools view -Oz -o {output_vcf_path}.vcf.gz {output_vcf_path}.bcf\n"
        script_content += f"tabix -p vcf {output_vcf_path}.vcf.gz\n"
        script_content += f"gunzip -c {output_vcf_path}.vcf.gz > {output_vcf_path}.vcf\n"
        script_content += f"rm -f {output_vcf_path}.bcf {output_vcf_path}.vcf.gz {output_vcf_path}.bcf.csi {output_vcf_path}.vcf.gz.tbi\n"
        script_content += "echo \"Variant calling and file processing completed.\"\n"
    
    script_path = '/data/Script_Site/tmp/vcf_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="vcf_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=vcf_script.sh"
    return response

@vcf_bp.route('/get_configurations_vcf', methods=['GET'])
@role_requis('superadmin')
def get_configurations_vcf():
    return jsonify(configurations_vcf)
    
@vcf_bp.route('/delete_config_vcf', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_vcf():
    index = request.json['index']
    try:
        configurations_vcf.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration not found.")
