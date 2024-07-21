from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import FullWorkflowConfiguration
from utils import role_requis

full_workflow_bp = Blueprint('full_workflow_bp', __name__)

configurations_full_workflow = []

@full_workflow_bp.route('/full_workflow', methods=['GET', 'POST'])
@role_requis('superadmin')
def full_workflow():
    if request.method == 'POST':
        base_output_dir = request.form['base_output_dir']
        input_dir = request.form['input_dir']
        ref_genome = request.form['ref_genome']
        qs_scores = request.form['qs_scores']
        cuda_device = request.form['cuda_device']
        model = request.form['model']
        kit_name = request.form['kit_name']
        vcf_output_file = request.form['vcf_output_file']
        
        configurations_full_workflow.append({
            "base_output_dir": base_output_dir,
            "input_dir": input_dir,
            "ref_genome": ref_genome,
            "qs_scores": qs_scores,
            "cuda_device": cuda_device,
            "model": model,
            "kit_name": kit_name,
            "vcf_output_file": vcf_output_file
        })

        new_config = FullWorkflowConfiguration(
            base_output_dir=base_output_dir,
            input_dir=input_dir,
            ref_genome=ref_genome,
            qs_scores=qs_scores,
            cuda_device=cuda_device,
            model=model,
            kit_name=kit_name,
            vcf_output_file=vcf_output_file
        )
        db.session.add(new_config)
        db.session.commit()
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('full_workflow.html')

@full_workflow_bp.route('/generate_full_workflow_script', methods=['GET'])
@role_requis('superadmin') 
def generate_full_workflow_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    
    for config in configurations_full_workflow:
        dorado_bin = "/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
        model_path = f"/home/grid/dorado-0.7.2-linux-x64/bin/{config['model']}"
        
        qs_scores_list = config['qs_scores'].split()
        for qscore in qs_scores_list:
            demultiplexed_dir = f"{config['base_output_dir']}/demultiplexed_q{qscore}"
            script_content += f"mkdir -p \"{demultiplexed_dir}\"\n"
            script_content += f"echo \"Starting Basecalling and Demultiplexing for Q-score {qscore}...\"\n"
            script_content += f"{dorado_bin} basecaller -x \"{config['cuda_device']}\" --min-qscore \"{qscore}\" --no-trim --emit-fastq {model_path} \"{config['input_dir']}\" | \\\n"
            script_content += f"{dorado_bin} demux --kit-name \"{config['kit_name']}\" --emit-fastq --output-dir \"{demultiplexed_dir}\"\n"
            script_content += "echo \"Processing complete for {config['input_dir']} with Q-score {qscore}\"\n"
            script_content += f"for fastq_file in \"{demultiplexed_dir}\"/*.fastq; do\n"
            script_content += f"    bam_file=\"${{fastq_file%.fastq}}.bam\"\n"
            script_content += f"    echo \"Aligning ${{fastq_file}} to reference genome...\"\n"
            script_content += f"    minimap2 -ax map-ont \"{config['ref_genome']}\" \"$fastq_file\" | samtools sort -o \"$bam_file\"\n"
            script_content += f"    samtools index \"$bam_file\"\n"
            script_content += f"    echo \"Alignment and BAM conversion completed for $bam_file\"\n"
            script_content += "done\n"

            vcf_filename = f"{config['vcf_output_file']}_q{qscore}.vcf"
            script_content += "echo \"Starting VCF generation...\"\n"
            script_content += f"samtools faidx \"{config['ref_genome']}\"\n"
            script_content += f"samtools index \"$bam_file\"\n"
            script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"$bam_file\" | bcftools call -mv -Ob -o \"{vcf_filename}.bcf\"\n"
            script_content += f"bcftools index \"{vcf_filename}.bcf\"\n"
            script_content += f"bcftools view -Oz -o \"{vcf_filename}.gz\" \"{vcf_filename}.bcf\"\n"
            script_content += f"tabix -p vcf \"{vcf_filename}.gz\"\n"
            script_content += f"gunzip -c \"{vcf_filename}.gz\" > \"{vcf_filename}\"\n"
            script_content += "echo \"VCF generation completed.\"\n"

    script_content += "echo \"All processes are complete.\"\n"
    return jsonify(script=script_content)

@full_workflow_bp.route('/download_full_workflow_script', methods=['GET'])
@role_requis('superadmin')
def download_full_workflow_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    
    for config in configurations_full_workflow:
        dorado_bin = "/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
        model_path = f"/home/grid/dorado-0.7.2-linux-x64/bin/{config['model']}"
        
        qs_scores_list = config['qs_scores'].split()
        for qscore in qs_scores_list:
            demultiplexed_dir = f"{config['base_output_dir']}/demultiplexed_q{qscore}"
            script_content += f"mkdir -p \"{demultiplexed_dir}\"\n"
            script_content += f"echo \"Starting Basecalling and Demultiplexing for Q-score {qscore}...\"\n"
            script_content += f"{dorado_bin} basecaller -x \"{config['cuda_device']}\" --min-qscore \"{qscore}\" --no-trim --emit-fastq {model_path} \"{config['input_dir']}\" | \\\n"
            script_content += f"{dorado_bin} demux --kit-name \"{config['kit_name']}\" --emit-fastq --output-dir \"{demultiplexed_dir}\"\n"
            script_content += "echo \"Processing complete for {config['input_dir']} with Q-score {qscore}\"\n"
            script_content += f"for fastq_file in \"{demultiplexed_dir}\"/*.fastq; do\n"
            script_content += f"    bam_file=\"${{fastq_file%.fastq}}.bam\"\n"
            script_content += f"    echo \"Aligning ${{fastq_file}} to reference genome...\"\n"
            script_content += f"    minimap2 -ax map-ont \"{config['ref_genome']}\" \"$fastq_file\" | samtools sort -o \"$bam_file\"\n"
            script_content += f"    samtools index \"$bam_file\"\n"
            script_content += f"    echo \"Alignment and BAM conversion completed for $bam_file\"\n"
            script_content += "done\n"

            vcf_filename = f"{config['vcf_output_file']}_q{qscore}.vcf"
            script_content += "echo \"Starting VCF generation...\"\n"
            script_content += f"samtools faidx \"{config['ref_genome']}\"\n"
            script_content += f"samtools index \"$bam_file\"\n"
            script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"$bam_file\" | bcftools call -mv -Ob -o \"{vcf_filename}.bcf\"\n"
            script_content += f"bcftools index \"{vcf_filename}.bcf\"\n"
            script_content += f"bcftools view -Oz -o \"{vcf_filename}.gz\" \"{vcf_filename}.bcf\"\n"
            script_content += f"tabix -p vcf \"{vcf_filename}.gz\"\n"
            script_content += f"gunzip -c \"{vcf_filename}.gz\" > \"{vcf_filename}\"\n"
            script_content += "echo \"VCF generation completed.\"\n"
            
    script_path = '/data/Script_Site/tmp/full_workflow_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
        
    response = make_response(send_file(script_path, as_attachment=True, download_name="full_workflow_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=full_workflow_script.sh"
    return response

@full_workflow_bp.route('/get_configurations_full_workflow', methods=['GET'])
@role_requis('superadmin')
def get_configurations_full_workflow():
    return jsonify(configurations_full_workflow)

@full_workflow_bp.route('/delete_config_full_workflow', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_full_workflow():
    index = request.json['index']
    try:
        db.session.delete(FullWorkflowConfiguration.query.get(index))
        db.session.commit()
        return jsonify(success=True, message="Configuration deleted successfully.")
    except Exception as e:
        return jsonify(success=False, message=str(e))
