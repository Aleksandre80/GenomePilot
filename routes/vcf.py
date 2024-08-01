from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationVCF, Workflow
import os
from utils import role_requis
import subprocess
import shlex
import json
from datetime import datetime
import re

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
        bam_filename = os.path.basename(config['bam_file'])
        bam_basename = bam_filename.replace('.bam', '')
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        
        # Extract q-score from BAM filename
        qscore = "X"  # Default to "X" if not found
        qscore_match = re.search(r'_q(\d+)', bam_basename)
        if qscore_match:
            qscore = qscore_match.group(1)
        
        vcf_directory = os.path.join(config['output_dir'], f"VCF_q{qscore}")
        log_file = f"{vcf_directory}/vcf_log.txt"
        report_file = f"{vcf_directory}/vcf_report.html"
        status_file = f"{vcf_directory}/vcf_status.txt"
        output_vcf_path = os.path.join(vcf_directory, f"{bam_basename}_{ref_basename}.vcf")

        script_content += f"mkdir -p \"{vcf_directory}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting VCF generation for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        script_content += f"samtools faidx \"{config['ref_genome']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"samtools index \"{config['bam_file']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished indexing {config['bam_file']}\" >> \"{log_file}\"\n"

        # Commandes pour générer le VCF sans logging des contenus intermédiaires
        script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"{config['bam_file']}\" | bcftools call -mv -Ob -o \"{output_vcf_path}.bcf\"\n"
        script_content += f"bcftools index \"{output_vcf_path}.bcf\"\n"
        script_content += f"bcftools view -Oz -o \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf\"\n"
        script_content += f"tabix -p vcf \"{output_vcf_path}.vcf.gz\"\n"
        script_content += f"gunzip -c \"{output_vcf_path}.vcf.gz\" > \"{output_vcf_path}\"\n"
        script_content += f"rm -f \"{output_vcf_path}.bcf\" \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf.csi\" \"{output_vcf_path}.vcf.gz.tbi\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished VCF processing for {config['bam_file']}\" >> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{output_vcf_path}\" ]; then\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>VCF Log Report</title></head><body><div class=\"log-container\"><h1>VCF Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"

    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)


@vcf_bp.route('/download_vcf_script', methods=['GET'])
@role_requis('superadmin')
def download_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        bam_filename = os.path.basename(config['bam_file'])
        bam_basename = bam_filename.replace('.bam', '')
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        
        # Extract q-score from BAM filename
        qscore = "X"  # Default to "X" if not found
        qscore_match = re.search(r'_q(\d+)', bam_basename)
        if qscore_match:
            qscore = qscore_match.group(1)
        
        vcf_directory = os.path.join(config['output_dir'], f"VCF_q{qscore}")
        log_file = f"{vcf_directory}/vcf_log.txt"
        report_file = f"{vcf_directory}/vcf_report.html"
        status_file = f"{vcf_directory}/vcf_status.txt"
        output_vcf_path = os.path.join(vcf_directory, f"{bam_basename}_{ref_basename}.vcf")

        script_content += f"mkdir -p \"{vcf_directory}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting VCF generation for BAM file {config['bam_file']}\" >> \"{log_file}\"\n"

        script_content += f"samtools faidx \"{config['ref_genome']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"samtools index \"{config['bam_file']}\" >> \"{log_file}\" 2>&1\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished indexing {config['bam_file']}\" >> \"{log_file}\"\n"

        # Commandes pour générer le VCF sans logging des contenus intermédiaires
        script_content += f"bcftools mpileup -Ou -f \"{config['ref_genome']}\" \"{config['bam_file']}\" | bcftools call -mv -Ob -o \"{output_vcf_path}.bcf\"\n"
        script_content += f"bcftools index \"{output_vcf_path}.bcf\"\n"
        script_content += f"bcftools view -Oz -o \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf\"\n"
        script_content += f"tabix -p vcf \"{output_vcf_path}.vcf.gz\"\n"
        script_content += f"gunzip -c \"{output_vcf_path}.vcf.gz\" > \"{output_vcf_path}\"\n"
        script_content += f"rm -f \"{output_vcf_path}.bcf\" \"{output_vcf_path}.vcf.gz\" \"{output_vcf_path}.bcf.csi\" \"{output_vcf_path}.vcf.gz.tbi\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished VCF processing for {config['bam_file']}\" >> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{output_vcf_path}\" ]; then\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>VCF Log Report</title></head><body><div class=\"log-container\"><h1>VCF Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"


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
def handle_vcf_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="VCF Generation", status="Running", start_time=datetime.utcnow(), output_dir=configurations_vcf[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/vcf_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_vcf[-1]['output_dir'] + "/vcf/vcf_status.txt"
            if os.path.exists(status_file):
                with open(status_file, 'r') as file:
                    status_info = file.read().strip()
                    status, end_time = status_info.split(' - ')
                    new_workflow.status = "Completed" if status == "completed" else "Failed"
                    new_workflow.end_time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
            else:
                new_workflow.status = "Failed"
                print("not found")
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

@vcf_bp.route('/history-vcf')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationVCF.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-vcf.html', configurations=configurations)