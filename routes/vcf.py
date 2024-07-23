from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationVCF, Workflow
import os
from utils import role_requis
import subprocess
import shlex

vcf_bp = Blueprint('vcf_bp', __name__)

configurations_vcf = []

@vcf_bp.route('/vcf_creator', methods=['GET', 'POST'])
@role_requis('superadmin') 
def vcf_creator():
    if request.method == 'POST':
        ref_genome_path = request.form['ref_genome']
        bam_file = request.form['bam_file']
        output_dir = request.form['output_dir']

        configurations_vcf.append({
            "ref_genome": ref_genome_path,
            "bam_file": bam_file,
            "output_dir": output_dir
        })
        
        configurations_vcf_db = ConfigurationVCF(
            ref_genome=ref_genome_path,
            bam_file=bam_file,
            output_dir=output_dir
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
        vcf_directory = os.path.join(config['output_dir'], "vcf")
        log_file = f"{vcf_directory}/vcf_log.txt"
        report_file = f"{vcf_directory}/vcf_report.html"

        script_content += f"mkdir -p \"{vcf_directory}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting VCF generation for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        output_vcf_path = os.path.join(vcf_directory, os.path.basename(config['output_vcf']))

        script_content += f"samtools faidx \"{config['ref_genome']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"samtools index \"{config['bam_file']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"{config['bam_file']}\" | bcftools call -mv -Ob -o \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools index \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools view -Oz -o \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"tabix -p vcf \"{output_vcf_path}.vcf.gz\" >> \"{log_file}\" 2>&1\n"
        script_content += f"gunzip -c \"{output_vcf_path}.vcf.gz\" > \"{output_vcf_path}.vcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"rm -f \"{output_vcf_path}.bcf\" \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf.csi\" \"{output_vcf_path}.vcf.gz.tbi\" >> \"{log_file}\" 2>&1\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Variant calling and file processing completed for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>VCF Log Report</title></head><body><div class=\"log-container\"><h1>VCF Log Report</h1>' > {report_file}\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> {report_file}\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> {report_file}\n"
    return jsonify(script=script_content)


@vcf_bp.route('/download_vcf_script', methods=['GET'])
@role_requis('superadmin')
def download_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        vcf_directory = os.path.join(config['output_dir'], "vcf")
        log_file = f"{vcf_directory}/vcf_log.txt"
        report_file = f"{vcf_directory}/vcf_report.html"

        script_content += f"mkdir -p \"{vcf_directory}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting VCF generation for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        output_vcf_path = os.path.join(vcf_directory, os.path.basename(config['output_vcf']))

        script_content += f"samtools faidx \"{config['ref_genome']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"samtools index \"{config['bam_file']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"{config['bam_file']}\" | bcftools call -mv -Ob -o \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools index \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"bcftools view -Oz -o \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"tabix -p vcf \"{output_vcf_path}.vcf.gz\" >> \"{log_file}\" 2>&1\n"
        script_content += f"gunzip -c \"{output_vcf_path}.vcf.gz\" > \"{output_vcf_path}.vcf\" >> \"{log_file}\" 2>&1\n"
        script_content += f"rm -f \"{output_vcf_path}.bcf\" \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf.csi\" \"{output_vcf_path}.vcf.gz.tbi\" >> \"{log_file}\" 2>&1\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Variant calling and file processing completed for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>VCF Log Report</title></head><body><div class=\"log-container\"><h1>VCF Log Report</h1>' > {report_file}\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> {report_file}\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> {report_file}\n"

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

@vcf_bp.route('/start_vcf_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="BAM Merge", status="Running")
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/bam_merge_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            # Assuming report_file path is stored in a way that it can be dynamically resolved
            report_file = configurations_vcf[-1]['output_dir'] + "/merge_report.html"
            if os.path.exists(report_file):
                new_workflow.status = "Completed"
            else:
                new_workflow.status = "Failed"
            
            db.session.commit()

        except Exception as e:
            print(f"Error: {e}")
            new_workflow.status = "Failed"
            db.session.commit()

        return jsonify(success=True, report=report_file)
    
    return jsonify(success=False, message="Invalid request method. Use POST.")