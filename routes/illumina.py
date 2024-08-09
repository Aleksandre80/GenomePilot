import subprocess
import shlex
import os
from flask import Blueprint, render_template, request, jsonify, make_response, send_file
from extensions import db
from models import ConfigurationIllumina, Workflow
from utils import role_requis
from datetime import datetime
import json

illumina_bp = Blueprint('illumina_bp', __name__)

configurations_illumina = []

@illumina_bp.route('/illumina', methods=['GET', 'POST'])
@role_requis('superadmin')
def illumina():
    if request.method == 'POST':
        input_dir = request.form['input_dir']
        output_dir = request.form['output_dir']
        ref_genome = request.form['ref_genome']
        
        if not all([input_dir, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_illumina.append({
            "input_dir": input_dir,
            "output_dir": output_dir,
            "ref_genome": ref_genome
        })
        
        configurations_illumina_db = ConfigurationIllumina(
            input_dir=input_dir,
            output_dir=output_dir,
            ref_genome=ref_genome
        )
        db.session.add(configurations_illumina_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('illumina.html')

@illumina_bp.route('/generate_illumina_script', methods=['GET'])
@role_requis('superadmin')
def generate_illumina_script():
    script_content = "#!/bin/bash\n\n"
    
    for config in configurations_illumina:
        output_dir = config['output_dir']
        log_file = f"{output_dir}/fastq_to_bam_log.txt"
        report_file = f"{output_dir}/fastq_to_bam_report.html"
        status_file = f"{output_dir}/fastq_to_bam_status.txt"
        reference_genome = config['ref_genome']
        input_dir = config['input_dir']

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting FASTQ to BAM conversion for directory {input_dir}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"

        # Scan for FASTQ pairs
        script_content += f"""
        for fastq_file in {input_dir}/*_R1_001.fastq; do
            base_name=$(basename "$fastq_file" "_R1_001.fastq")
            fastq_r1="{input_dir}/${{base_name}}_R1_001.fastq"
            fastq_r2="{input_dir}/${{base_name}}_R2_001.fastq"
            if [[ -f "$fastq_r2" ]]; then
                echo "Found pair: $fastq_r1 and $fastq_r2" >> "{log_file}"
                
                output_bam="{output_dir}/${{base_name}}_aligned.bam"
                output_sorted_bam="{output_dir}/${{base_name}}_sorted.bam"
                
                # Minimap2 alignment
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Aligning $fastq_r1 and $fastq_r2 to reference genome..." >> "{log_file}"
                minimap2 -ax sr "{reference_genome}" "$fastq_r1" "$fastq_r2" > "${{output_bam}}" 2>> "{log_file}"
                
                # Convert SAM to BAM and sort BAM
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Converting and sorting BAM file..." >> "{log_file}"
                samtools view -Sb "${{output_bam}}" | samtools sort -o "${{output_sorted_bam}}" 2>> "{log_file}"
                
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Processed $base_name" >> "{log_file}"
            else
                echo "Error: Could not find corresponding R2 file for $fastq_r1" >> "{log_file}"
            fi
        done

        # Merging BAM files for each patient
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Merging sorted BAM files..." >> "{log_file}"
        for patient_id in $(ls {output_dir}/*_sorted.bam | sed -n 's/.*_\\(S[0-9]\\+\\)_L.*/\\1/p' | sort | uniq); do
            bam_files=$(ls {output_dir}/*_${{patient_id}}_*_sorted.bam)
            if [[ -n "$bam_files" ]]; then
                merged_bam="{output_dir}/${{patient_id}}_merged.bam"
                echo "Merging BAM files for $patient_id: $bam_files" >> "{log_file}"
                samtools merge "$merged_bam" $bam_files 2>> "{log_file}"
                
                # Index merged BAM
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Indexing merged BAM file..." >> "{log_file}"
                samtools index "$merged_bam" 2>> "{log_file}"
                
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Merging and indexing completed for patient $patient_id" >> "{log_file}"
            else
                echo "No BAM files found for $patient_id" >> "{log_file}"
            fi
        done
        """

        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - All processes are complete.\" >> \"{log_file}\"\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>FASTQ to BAM Log Report</title></head><body><div class=\"log-container\"><h1>FASTQ to BAM Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)
    
    return jsonify(script=escaped_script_content)




@illumina_bp.route('/download_illumina_script', methods=['GET'])
@role_requis('superadmin')
def download_coberage_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_illumina:
        output_dir = config['output_dir']
        log_file = f"{output_dir}/fastq_to_bam_log.txt"
        report_file = f"{output_dir}/fastq_to_bam_report.html"
        status_file = f"{output_dir}/fastq_to_bam_status.txt"
        reference_genome = config['ref_genome']
        input_dir = config['input_dir']

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting FASTQ to BAM conversion for directory {input_dir}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"

        # Scan for FASTQ pairs
        script_content += f"""
        for fastq_file in {input_dir}/*_R1_001.fastq; do
            base_name=$(basename "$fastq_file" "_R1_001.fastq")
            fastq_r1="{input_dir}/${{base_name}}_R1_001.fastq"
            fastq_r2="{input_dir}/${{base_name}}_R2_001.fastq"
            if [[ -f "$fastq_r2" ]]; then
                echo "Found pair: $fastq_r1 and $fastq_r2" >> "{log_file}"
                
                output_bam="{output_dir}/${{base_name}}_aligned.bam"
                output_sorted_bam="{output_dir}/${{base_name}}_sorted.bam"
                
                # Minimap2 alignment
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Aligning $fastq_r1 and $fastq_r2 to reference genome..." >> "{log_file}"
                minimap2 -ax sr "{reference_genome}" "$fastq_r1" "$fastq_r2" > "${{output_bam}}" 2>> "{log_file}"
                
                # Convert SAM to BAM and sort BAM
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Converting and sorting BAM file..." >> "{log_file}"
                samtools view -Sb "${{output_bam}}" | samtools sort -o "${{output_sorted_bam}}" 2>> "{log_file}"
                
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Processed $base_name" >> "{log_file}"
            else
                echo "Error: Could not find corresponding R2 file for $fastq_r1" >> "{log_file}"
            fi
        done

        # Merging BAM files for each patient
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Merging sorted BAM files..." >> "{log_file}"
        for patient_id in $(ls {output_dir}/*_sorted.bam | sed -n 's/.*_\\(S[0-9]\\+\\)_L.*/\\1/p' | sort | uniq); do
            bam_files=$(ls {output_dir}/*_${{patient_id}}_*_sorted.bam)
            if [[ -n "$bam_files" ]]; then
                merged_bam="{output_dir}/${{patient_id}}_merged.bam"
                echo "Merging BAM files for $patient_id: $bam_files" >> "{log_file}"
                samtools merge "$merged_bam" $bam_files 2>> "{log_file}"
                
                # Index merged BAM
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Indexing merged BAM file..." >> "{log_file}"
                samtools index "$merged_bam" 2>> "{log_file}"
                
                echo "$(date '+%Y-%m-%d %H:%M:%S') - Merging and indexing completed for patient $patient_id" >> "{log_file}"
            else
                echo "No BAM files found for $patient_id" >> "{log_file}"
            fi
        done
        """

        script_content += f"if [ $? -eq 0 ]; then\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - All processes are complete.\" >> \"{log_file}\"\n"

        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>FASTQ to BAM Log Report</title></head><body><div class=\"log-container\"><h1>FASTQ to BAM Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"


    
    script_path = '/data/Script_Site/tmp/illumina_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="illumina_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=illumina_script.sh"
    return response



@illumina_bp.route('/get_configurations_illumina', methods=['GET'])
@role_requis('superadmin')
def get_configurations_illumina():
    return jsonify(configurations_illumina)

@illumina_bp.route('/delete_config_illumina', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_illumina():
    index = request.json['index']
    try:
        configurations_illumina.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration does not exist.")

@illumina_bp.route('/start_illumina_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_illumina_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="illumina", status="Running", start_time=datetime.utcnow(), output_dir=configurations_illumina[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/illumina_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = os.path.join(configurations_illumina[-1]['output_dir'], "illumina", "illumina_status.txt")
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


@illumina_bp.route('/history-illumina')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationIllumina.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-illumina.html', configurations=configurations)