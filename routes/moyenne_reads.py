from flask import Blueprint, render_template, request, jsonify, redirect, url_for, make_response, send_file
from extensions import db
from models import ConfigurationMoyenneReads, Workflow
from utils import role_requis
import json
import os
import subprocess
import shlex
from datetime import datetime

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
def generate_moyenne_reads_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_moyenne_reads:
        bam_filename = os.path.basename(config['bam_file'])
        bam_basename = bam_filename.replace('.bam', '')
        bed_filename = os.path.basename(config['bed_file']).replace('.bed', '')

        output_dir = os.path.join(config['output_dir'], "Moyenne-Reads")
        log_file = f"{output_dir}/moyenne_reads_log.txt"
        report_file = f"{output_dir}/moyenne_reads_report.html"
        status_file = f"{output_dir}/moyenne_reads_status.txt"
        coverage_txt = f"{output_dir}/{bam_basename}_{bed_filename}.txt"
        coverage_csv = f"{output_dir}/{bam_basename}_{bed_filename}.csv"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting Calculating reads average for BAM file {config['bam_file']} and BED file {config['bed_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour calculer la moyenne des reads
        script_content += f"bedtools coverage -a \"{config['bed_file']}\" -b \"{config['bam_file']}\" > \"{coverage_txt}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{coverage_txt}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Calculating reads average completed and {bam_basename}_{bed_filename}.txt generated.\" >> \"{log_file}\"\n"
        
        # Conversion de coverage.txt en coverage.csv
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Converting {bam_basename}_{bed_filename}.txt to {bam_basename}_{bed_filename}.csv...\" >> \"{log_file}\"\n"
        script_content += f"    awk 'BEGIN {{ OFS=\",\"; print \"chrom\",\"start\",\"end\",\"num_reads\",\"bases_covered\",\"coverage\" }} {{ print $1,$2,$3,$4,$5,$6 }}' \"{coverage_txt}\" > \"{coverage_csv}\" 2>> \"{log_file}\"\n"
        
        script_content += f"    if [ -f \"{coverage_csv}\" ]; then\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Conversion to {bam_basename}_{bed_filename}.csv completed.\" >> \"{log_file}\"\n"
        script_content += f"        echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"    else\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Conversion to {bam_basename}_{bed_filename}.csv failed.\" >> \"{log_file}\"\n"
        script_content += f"        echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"    fi\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Calculating reads average failed. {bam_basename}_{bed_filename}.txt not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Reads Average Log Report</title></head><body><div class=\"log-container\"><h1>Reads Average Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)





@moyenne_reads_bp.route('/download_moyenne_reads_script', methods=['GET'])
@role_requis('superadmin')
def download_moyenne_reads_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_moyenne_reads:
        bam_filename = os.path.basename(config['bam_file'])
        bam_basename = bam_filename.replace('.bam', '')
        bed_filename = os.path.basename(config['bed_file']).replace('.bed', '')

        output_dir = os.path.join(config['output_dir'], "Moyenne-Reads")
        log_file = f"{output_dir}/moyenne_reads_log.txt"
        report_file = f"{output_dir}/moyenne_reads_report.html"
        status_file = f"{output_dir}/moyenne_reads_status.txt"
        coverage_txt = f"{output_dir}/{bam_basename}_{bed_filename}.txt"
        coverage_csv = f"{output_dir}/{bam_basename}_{bed_filename}.csv"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting Calculating reads average for BAM file {config['bam_file']} and BED file {config['bed_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"
        
        # Commande pour calculer la moyenne des reads
        script_content += f"bedtools coverage -a \"{config['bed_file']}\" -b \"{config['bam_file']}\" > \"{coverage_txt}\" 2>> \"{log_file}\"\n"
        
        script_content += f"if [ -f \"{coverage_txt}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Calculating reads average completed and {bam_basename}_{bed_filename}.txt generated.\" >> \"{log_file}\"\n"
        
        # Conversion de coverage.txt en coverage.csv
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Converting {bam_basename}_{bed_filename}.txt to {bam_basename}_{bed_filename}.csv...\" >> \"{log_file}\"\n"
        script_content += f"    awk 'BEGIN {{ OFS=\",\"; print \"chrom\",\"start\",\"end\",\"num_reads\",\"bases_covered\",\"coverage\" }} {{ print $1,$2,$3,$4,$5,$6 }}' \"{coverage_txt}\" > \"{coverage_csv}\" 2>> \"{log_file}\"\n"
        
        script_content += f"    if [ -f \"{coverage_csv}\" ]; then\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Conversion to {bam_basename}_{bed_filename}.csv completed.\" >> \"{log_file}\"\n"
        script_content += f"        echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"    else\n"
        script_content += f"        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Conversion to {bam_basename}_{bed_filename}.csv failed.\" >> \"{log_file}\"\n"
        script_content += f"        echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"    fi\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Calculating reads average failed. {bam_basename}_{bed_filename}.txt not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Reads Average Log Report</title></head><body><div class=\"log-container\"><h1>Reads Average Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class='log-entry'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
        
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

@moyenne_reads_bp.route('/start_moyenne_reads_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_moyenne_reads_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Calculating Reads Average", status="Running", start_time=datetime.utcnow(), output_dir=configurations_moyenne_reads[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/moyenne_reads_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = configurations_moyenne_reads[-1]['output_dir'] + "/moyenne_reads_status.txt"
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
