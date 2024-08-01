from flask import Blueprint, render_template, request, jsonify, flash, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationBasecalling, Workflow
from utils import role_requis
import json
import os
import subprocess
import shlex
from datetime import datetime

basecalling_bp = Blueprint('basecalling_bp', __name__)

configurations_basecalling = []

@basecalling_bp.route('/basecalling', methods=['GET', 'POST'])
@role_requis('superadmin') 
def basecalling():
    if request.method == 'POST':
        base_output_dir = request.form['base_output_dir']
        input_dir = request.form['input_dir']
        ref_genome_path = request.form['ref_genome']
        qs_scores = request.form['qs_scores']
        cuda_device = request.form['cuda_device']
        model = request.form['model']
        kit_name = request.form['kit_name']

        if not all([base_output_dir, input_dir, ref_genome_path, qs_scores, cuda_device, model, kit_name]):
            flash('Please fill all fields before adding a configuration.', 'error')
            return redirect(url_for('basecalling_bp.basecalling'))
        
        configurations_basecalling.append({
            "base_output_dir": base_output_dir,
            "input_dir": input_dir,
            "ref_genome": ref_genome_path,
            "qs_scores": qs_scores,
            "cuda_device": cuda_device,
            "model": model,
            "kit_name": kit_name
        })
        
        configurations_basecalling_db = ConfigurationBasecalling(
            base_output_dir=base_output_dir,
            input_dir=input_dir,
            ref_genome=ref_genome_path,
            qs_scores=qs_scores,
            cuda_device=cuda_device,
            model=model,
            kit_name=kit_name
        )
        db.session.add(configurations_basecalling_db)
        db.session.commit()
        
        flash('Configuration added successfully.', 'success')
        return jsonify({'success': True, 'configurations': configurations_basecalling}), 200

    return render_template('index.html', configurations=configurations_basecalling)

@basecalling_bp.route('/generate_script', methods=['GET'])
@role_requis('superadmin')
def generate_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_basecalling:
        log_file = f"{config['base_output_dir']}/basecalling_log.txt"
        report_file = f"{config['base_output_dir']}/basecalling_report.html"
        status_file = f"{config['base_output_dir']}/basecalling_status.txt"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting basecalling for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        
        qs_scores_list = config['qs_scores'].split()
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        model_basename = os.path.basename(config['model']).replace('.bin', '')
        kit_name = config['kit_name']
        barcode_name = config.get('barcode_name', 'default_barcode')  # Assuming barcode_name is provided in config

        for qscore in qs_scores_list:
            output_dir = f"${{BASE_OUTPUT_DIR}}/demultiplexed_q{qscore}"
            script_content += f"BASE_OUTPUT_DIR=\"{config['base_output_dir']}\"\n"
            script_content += "mkdir -p \"${BASE_OUTPUT_DIR}\"\n"
            script_content += f"""
DORADO_BIN="/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
MODEL_PATH="/home/grid/dorado-0.7.2-linux-x64/bin/{config['model']}"
REF_GENOME="{config['ref_genome']}"
INPUT_DIR="{config['input_dir']}"
OUTPUT_DIR="{output_dir}"
mkdir -p "${{OUTPUT_DIR}}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting basecalling and demultiplexing for Q-score {qscore}" >> \"{log_file}\"
${{DORADO_BIN}} basecaller -x "{config['cuda_device']}" --min-qscore "{qscore}" --no-trim --emit-fastq ${{MODEL_PATH}} ${{INPUT_DIR}} | \\
${{DORADO_BIN}} demux --kit-name "{kit_name}" --emit-fastq --output-dir "${{OUTPUT_DIR}}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing complete for {config['input_dir']} with Q-score {qscore}" >> \"{log_file}\"
"""
            script_content += f"for fastq_file in \"${{OUTPUT_DIR}}\"/*.fastq; do\n"
            script_content += f"    bam_file=\"{ref_basename}_{model_basename}_{kit_name}_{barcode_name}_q{qscore}.bam\"\n"
            script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Aligning ${{fastq_file}} to reference genome...\" >> \"{log_file}\"\n"
            script_content += f"    minimap2 -ax map-ont \"{config['ref_genome']}\" \"$fastq_file\" | samtools sort -o \"${{OUTPUT_DIR}}/$bam_file\"\n"
            script_content += f"    samtools index \"${{OUTPUT_DIR}}/$bam_file\"\n"
            script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Alignment and BAM conversion completed for ${{bam_file}}\" >> \"{log_file}\"\n"
            script_content += "done\n"
    
    script_content += f"if [ $? -eq 0 ]; then\n"
    script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
    script_content += f"else\n"
    script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
    script_content += f"fi\n"
    
    script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - All processes are complete.\" >> \"{log_file}\"\n"
    
    # Generate HTML report
    script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Basecalling Log Report</title></head><body><div class=\"log-container\"><h1>Basecalling Log Report</h1>' > \"{report_file}\"\n"
    script_content += f"while IFS= read -r line; do\n"
    script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
    script_content += f"done < \"{log_file}\"\n"
    script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)





@basecalling_bp.route('/download_basecalling_script', methods=['GET'])
@role_requis('superadmin')
def download_basecalling_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_basecalling:
        log_file = f"{config['base_output_dir']}/basecalling_log.txt"
        report_file = f"{config['base_output_dir']}/basecalling_report.html"
        status_file = f"{config['base_output_dir']}/basecalling_status.txt"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting basecalling for input directory {config['input_dir']}\" >> \"{log_file}\"\n"
        
        qs_scores_list = config['qs_scores'].split()
        ref_basename = os.path.basename(config['ref_genome']).replace('.fa', '')
        model_basename = os.path.basename(config['model']).replace('.bin', '')
        kit_name = config['kit_name']
        barcode_name = config.get('barcode_name', 'default_barcode')  # Assuming barcode_name is provided in config

        for qscore in qs_scores_list:
            output_dir = f"${{BASE_OUTPUT_DIR}}/demultiplexed_q{qscore}"
            script_content += f"BASE_OUTPUT_DIR=\"{config['base_output_dir']}\"\n"
            script_content += "mkdir -p \"${BASE_OUTPUT_DIR}\"\n"
            script_content += f"""
DORADO_BIN="/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
MODEL_PATH="/home/grid/dorado-0.7.2-linux-x64/bin/{config['model']}"
REF_GENOME="{config['ref_genome']}"
INPUT_DIR="{config['input_dir']}"
OUTPUT_DIR="{output_dir}"
mkdir -p "${{OUTPUT_DIR}}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting basecalling and demultiplexing for Q-score {qscore}" >> \"{log_file}\"
${{DORADO_BIN}} basecaller -x "{config['cuda_device']}" --min-qscore "{qscore}" --no-trim --emit-fastq ${{MODEL_PATH}} ${{INPUT_DIR}} | \\
${{DORADO_BIN}} demux --kit-name "{kit_name}" --emit-fastq --output-dir "${{OUTPUT_DIR}}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing complete for {config['input_dir']} with Q-score {qscore}" >> \"{log_file}\"
"""
            script_content += f"for fastq_file in \"${{OUTPUT_DIR}}\"/*.fastq; do\n"
            script_content += f"    bam_file=\"{ref_basename}_{model_basename}_{kit_name}_{barcode_name}_q{qscore}.bam\"\n"
            script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Aligning ${{fastq_file}} to reference genome...\" >> \"{log_file}\"\n"
            script_content += f"    minimap2 -ax map-ont \"{config['ref_genome']}\" \"$fastq_file\" | samtools sort -o \"${{OUTPUT_DIR}}/$bam_file\"\n"
            script_content += f"    samtools index \"${{OUTPUT_DIR}}/$bam_file\"\n"
            script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Alignment and BAM conversion completed for ${{bam_file}}\" >> \"{log_file}\"\n"
            script_content += "done\n"
    
    script_content += f"if [ $? -eq 0 ]; then\n"
    script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
    script_content += f"else\n"
    script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
    script_content += f"fi\n"
    
    script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - All processes are complete.\" >> \"{log_file}\"\n"
    
    # Generate HTML report
    script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Basecalling Log Report</title></head><body><div class=\"log-container\"><h1>Basecalling Log Report</h1>' > \"{report_file}\"\n"
    script_content += f"while IFS= read -r line; do\n"
    script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
    script_content += f"done < \"{log_file}\"\n"
    script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    script_path = '/data/Script_Site/tmp/basecalling_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="basecalling_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=basecalling_script.sh"
    return response

@basecalling_bp.route('/get_configurations_basecalling', methods=['GET'])
@role_requis('superadmin')
def get_configurations_basecalling():
    return jsonify(configurations_basecalling)

@basecalling_bp.route('/delete_config_basecalling', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_basecalling():
    index = request.json['index']
    try:
        configurations_basecalling.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration not found.")
    
@basecalling_bp.route('/start_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_basecalling_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Basecalling", status="Running", start_time=datetime.utcnow(), output_dir=configurations_basecalling[-1]['base_output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/basecalling_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_basecalling[-1]['base_output_dir'] + "/basecalling_status.txt"
            if os.path.exists(status_file):
                with open(status_file, 'r') as file:
                    status_info = file.read().strip()
                    status, end_time = status_info.split(' - ')
                    new_workflow.status = "Completed" if status == "completed" else "Failed"
                    new_workflow.end_time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
            else:
                new_workflow.status = "Failed"
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

@basecalling_bp.route('/history-basecalling')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationBasecalling.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-basecalling.html', configurations=configurations)