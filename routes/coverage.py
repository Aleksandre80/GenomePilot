import subprocess
import shlex
import os
from flask import Blueprint, render_template, request, jsonify, make_response, send_file
from extensions import db
from models import ConfigurationCoverage, Workflow
from utils import role_requis
from datetime import datetime
import json

coverage_bp = Blueprint('coverage_bp', __name__)

configurations_coverage = []

@coverage_bp.route('/coverage', methods=['GET', 'POST'])
@role_requis('superadmin')
def coverage():
    if request.method == 'POST':
        input_file = request.form['input_file']
        output_dir = request.form['output_dir']
        bed_file = request.form['bed_file']
        
        if not all([input_file, output_dir]):
            return jsonify(success=False, message="Please specify both input and output directories.")
        
        configurations_coverage.append({
            "input_file": input_file,
            "output_dir": output_dir,
            "bed_file": bed_file
        })
        
        configurations_coverage_db = ConfigurationCoverage(
            input_file=input_file,
            output_dir=output_dir,
            bed_file=bed_file
        )
        db.session.add(configurations_coverage_db)
        db.session.commit()
        
        return jsonify(success=True, message="Configuration added successfully.")
    return render_template('coverage.html')

@coverage_bp.route('/generate_coverage_script', methods=['GET'])
@role_requis('superadmin')
def generate_coverage_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_coverage:
        bam_filename = os.path.basename(config['input_file'])
        bam_basename = bam_filename.replace('.bam', '')
        bed_filename = os.path.basename(config['bed_file']).replace('.bed', '')

        output_dir = os.path.join(config['output_dir'], "Coverage_Calculation")
        log_file = f"{output_dir}/coverage_log.txt"
        report_file = f"{output_dir}/coverage_report.html"
        status_file = f"{output_dir}/coverage_status.txt"
        coverage_output = f"{output_dir}/{bam_basename}_{bed_filename}_coverage_output.txt"
        avg_coverage_output = f"{output_dir}/{bam_basename}_{bed_filename}_average_coverages.txt"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting coverage calculation for BAM file {config['input_file']} and BED file {config['bed_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"

        # Commande pour calculer la couverture
        script_content += f"bedtools coverage -d -a \"{config['bed_file']}\" -b \"{config['input_file']}\" > \"{coverage_output}\" 2>> \"{log_file}\"\n"
        
        # Commande pour calculer la moyenne de couverture
        script_content += f"awk 'BEGIN {{ FS = \"\\t\"; OFS = \"\\t\"; }} " + \
                          f"{{ region = $1 \":\" $2 \"-\" $3 \" \" $4; coverage[region] += $6; count[region] += 1; }} " + \
                          f"END {{ for (r in coverage) {{ avg_coverage = coverage[r] / count[r]; print r, avg_coverage; }} }}' " + \
                          f"\"{coverage_output}\" > \"{avg_coverage_output}\" 2>> \"{log_file}\"\n"

        script_content += f"if [ -f \"{avg_coverage_output}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Coverage calculation completed and {avg_coverage_output} generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Coverage calculation failed. {avg_coverage_output} not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Coverage Calculation Log Report</title></head><body><div class=\"log-container\"><h1>Coverage Calculation Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class=\'log-entry\'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    # Échapper les caractères spéciaux pour JSON
    escaped_script_content = json.dumps(script_content)

    return jsonify(script=escaped_script_content)


@coverage_bp.route('/download_coverage_script', methods=['GET'])
@role_requis('superadmin')
def download_coberage_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_coverage:
        bam_filename = os.path.basename(config['input_file'])
        bam_basename = bam_filename.replace('.bam', '')
        bed_filename = os.path.basename(config['bed_file']).replace('.bed', '')

        output_dir = os.path.join(config['output_dir'], "Coverage_Calculation")
        log_file = f"{output_dir}/coverage_log.txt"
        report_file = f"{output_dir}/coverage_report.html"
        status_file = f"{output_dir}/coverage_status.txt"
        coverage_output = f"{output_dir}/{bam_basename}_{bed_filename}_coverage_output.txt"
        avg_coverage_output = f"{output_dir}/{bam_basename}_{bed_filename}_average_coverages.txt"

        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting coverage calculation for BAM file {config['input_file']} and BED file {config['bed_file']}\" >> \"{log_file}\"\n"
        script_content += f"mkdir -p \"{output_dir}\"\n"
        script_content += f"echo \"$(date '+%Y-%m-%d %H:%M:%S') - Output directory created.\" >> \"{log_file}\"\n"

        # Commande pour calculer la couverture
        script_content += f"bedtools coverage -d -a \"{config['bed_file']}\" -b \"{config['input_file']}\" > \"{coverage_output}\" 2>> \"{log_file}\"\n"
        
        # Commande pour calculer la moyenne de couverture
        script_content += f"awk 'BEGIN {{ FS = \"\\t\"; OFS = \"\\t\"; }} " + \
                          f"{{ region = $1 \":\" $2 \"-\" $3 \" \" $4; coverage[region] += $6; count[region] += 1; }} " + \
                          f"END {{ for (r in coverage) {{ avg_coverage = coverage[r] / count[r]; print r, avg_coverage; }} }}' " + \
                          f"\"{coverage_output}\" > \"{avg_coverage_output}\" 2>> \"{log_file}\"\n"

        script_content += f"if [ -f \"{avg_coverage_output}\" ]; then\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Coverage calculation completed and {avg_coverage_output} generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"completed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"else\n"
        script_content += f"    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Coverage calculation failed. {avg_coverage_output} not generated.\" >> \"{log_file}\"\n"
        script_content += f"    echo \"failed - $(date '+%Y-%m-%d %H:%M:%S')\" > \"{status_file}\"\n"
        script_content += f"fi\n"
        
        # Generate HTML report
        script_content += f"echo '<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>Coverage Calculation Log Report</title></head><body><div class=\"log-container\"><h1>Coverage Calculation Log Report</h1>' > \"{report_file}\"\n"
        script_content += f"while IFS= read -r line; do\n"
        script_content += f"    echo \"<div class=\'log-entry\'>\"$line\"</div>\" >> \"{report_file}\"\n"
        script_content += f"done < \"{log_file}\"\n"
        script_content += f"echo '</div></body></html>' >> \"{report_file}\"\n"
    
    script_path = '/data/Script_Site/tmp/coverage_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    response = make_response(send_file(script_path, as_attachment=True, download_name="coverage_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=coverage_script.sh"
    return response



@coverage_bp.route('/get_configurations_coverage', methods=['GET'])
@role_requis('superadmin')
def get_configurations_coverage():
    return jsonify(configurations_coverage)

@coverage_bp.route('/delete_config_coverage', methods=['POST'])
@role_requis('superadmin')
def delete_configuration_coverage():
    index = request.json['index']
    try:
        configurations_coverage.pop(index)
        return jsonify(success=True, message="Configuration deleted successfully.")
    except IndexError:
        return jsonify(success=False, message="Configuration does not exist.")

@coverage_bp.route('/start_coverage_script', methods=['GET', 'POST'])
@role_requis('superadmin')
def handle_coverage_script():
    if request.method == 'POST':
        new_workflow = Workflow(name="Coverage Calculation", status="Running", start_time=datetime.utcnow(), output_dir=configurations_coverage[-1]['output_dir'])
        db.session.add(new_workflow)
        db.session.commit()
        
        try:
            script_path = '/data/Script_Site/tmp/coverage_script.sh'
            script_command = f"bash {script_path}"
            
            process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()
            
            status_file = os.path.join(configurations_coverage[-1]['output_dir'], "Coverage_Calculation", "coverage_status.txt")
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


@coverage_bp.route('/history-coverage')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationCoverage.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-coverage.html', configurations=configurations)