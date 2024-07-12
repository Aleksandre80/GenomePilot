import os
from flask import Flask, render_template, request, redirect, url_for, flash, jsonify, session, send_file, make_response, copy_current_request_context
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from functools import wraps
from flask_migrate import Migrate
import subprocess
from flask_socketio import SocketIO, emit
import threading
import shlex

app = Flask(__name__)
app.secret_key = 'votre_clé_secrète_ici'
app.config['UPLOAD_FOLDER'] = 'path_to_save'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///configurations.db'
db = SQLAlchemy(app)
migrate = Migrate(app, db)
socketio = SocketIO(app, cors_allowed_origins="*")

configurations_basecalling = []
configurations_merge = []
configurations_vcf = []
configurations_full_workflow = []

class ConfigurationBasecalling(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    base_output_dir = db.Column(db.String(120), nullable=False)
    input_dir = db.Column(db.String(120), nullable=False)
    ref_genome = db.Column(db.String(120), nullable=False)
    qs_scores = db.Column(db.String(120), nullable=False)
    cuda_device = db.Column(db.String(120), nullable=False)
    model = db.Column(db.String(120), nullable=False)
    kit_name = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationMerge(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_dir = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationVCF(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ref_genome = db.Column(db.String(120), nullable=False)
    bam_file = db.Column(db.String(120), nullable=False)
    output_vcf = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class FullWorkflowConfiguration(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    base_output_dir = db.Column(db.String(120), nullable=False)
    input_dir = db.Column(db.String(120), nullable=False)
    ref_genome = db.Column(db.String(120), nullable=False)
    qs_scores = db.Column(db.String(120), nullable=False)
    cuda_device = db.Column(db.String(120), nullable=False)
    model = db.Column(db.String(120), nullable=False)
    kit_name = db.Column(db.String(120), nullable=False)
    vcf_output_file = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)

class Workflow(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False)
    launch_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    status = db.Column(db.String(20), nullable=False, default='Pending')

    def __repr__(self):
        return f"Workflow('{self.name}', '{self.launch_date}', '{self.status}')"

def role_requis(*roles_requis):
    def wrapper(fn):
        @wraps(fn)
        def decorated_view(*args, **kwargs):
            if get_role_utilisateur() not in roles_requis:
                flash("Accès refusé. Veuillez vous connecter avec un compte autorisé.", "warning")
                return redirect(url_for('login'))
            return fn(*args, **kwargs)
        return decorated_view
    return wrapper

@app.route('/', methods=['GET', 'POST'])
@role_requis('superadmin') 
def accueil():
    return render_template('workflow.html')

@app.route('/status', methods=['GET', 'POST'])
@role_requis('superadmin')
def status():
    workflows = Workflow.query.order_by(Workflow.launch_date.desc()).all()
    return render_template('status.html', workflows=workflows)

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        password = request.form['password']
        role = "user"
        if password == "Aleksandre":
            role = "superadmin"
        session['password'] = password
        session['role'] = role
        
        if role == "admin" or role == "superadmin":
            return redirect(url_for('accueil'))
        else:
            return redirect(url_for('login'))
            
    return render_template('login.html')

def get_role_utilisateur():
    return session.get('role', 'user')
    
@app.route('/full_workflow', methods=['GET', 'POST'])
@role_requis('superadmin')
def full_workflow():
    if request.method == 'POST':
        
        base_output_dir=request.form['base_output_dir']
        input_dir=request.form['input_dir']
        ref_genome=request.form['ref_genome']
        qs_scores=request.form['qs_scores']
        cuda_device=request.form['cuda_device']
        model=request.form['model']
        kit_name=request.form['kit_name']
        vcf_output_file=request.form['vcf_output_file']
        
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

        print(f"full_workflow", configurations_full_workflow)
        
        
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


@app.route('/generate_full_workflow_script', methods=['GET'])
@role_requis('superadmin') 
def generate_full_workflow_script():
    if not configurations_full_workflow:  # Vérifie si la liste est vide
        return jsonify(success=False, message="No configurations available")

    script_content = "#!/bin/bash\n\nsource ~/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    
    for config in configurations_full_workflow:
        # Définir les chemins des outils et des dossiers pour chaque qscore
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

            # Générer le VCF pour chaque qscore
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

@app.route('/vcf_creator', methods=['GET', 'POST'])
@role_requis('superadmin') 
def vcf_creator():
    if request.method == 'POST':
        # Extract file name and save the file
        ref_genome_path = request.form['ref_genome']
        bam_file = request.form['bam_file']
        output_vcf = request.form['output_vcf']

        # Append the path instead of the FileStorage object
        configurations_vcf.append({
            "ref_genome": ref_genome_path,  # Store path only
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


@app.route('/generate_vcf_script', methods=['GET'])
@role_requis('superadmin') 
def generate_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        # Création d'un dossier VCF si non existant
        output_directory = os.path.dirname(config['output_vcf'])
        vcf_directory = f"{output_directory}/vcf"
        script_content += f"mkdir -p {vcf_directory}\n"  # S'assure que le dossier vcf existe
        
        output_vcf_path = f"{vcf_directory}/{os.path.basename(config['output_vcf'])}"  # Chemin modifié pour placer les fichiers dans le dossier vcf
        
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

@app.route('/download_vcf_script', methods=['GET'])
@role_requis('superadmin')
def download_vcf_script():
    script_content = "#!/bin/bash\n\nsource /home/grid/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        # Création d'un dossier VCF dans le même répertoire que le fichier BAM
        bam_directory = os.path.dirname(config['bam_file'])
        vcf_directory = os.path.join(bam_directory, "vcf")
        script_content += f"mkdir -p {vcf_directory}\n"  # S'assure que le dossier vcf existe

        # Chemin modifié pour placer les fichiers VCF dans le dossier vcf
        base_vcf_filename = os.path.basename(config['output_vcf'])
        output_vcf_path = os.path.join(vcf_directory, base_vcf_filename)

        # Processus de création des fichiers VCF
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


@socketio.on('start_vcf_script', namespace='/vcf')
def handle_vcf_script():
    print("Received start script event from client.")
    emit("yoooo", namespace='/vcf')
    new_workflow = Workflow(name="VCF Creator", status="Running")
    db.session.add(new_workflow)
    db.session.commit()
    script_command = "bash /data/Script_Site/tmp/vcf_script.sh"
    try:
        process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        if stderr:
            emit('script_output', {'message': stderr, 'type': 'stderr'}, namespace='/vcf')
            new_workflow.status = "Failed"
            emit('script_error', {'status': "Failed"}, namespace='/vcf')
        if stdout:
            for line in stdout.splitlines():
                emit('script_output', {'message': line, 'type': 'stdout'}, namespace='/vcf')
                new_workflow.status = "Completed"
                emit('script_error', {'status': "Completed"}, namespace='/vcf')


        db.session.commit()

    except Exception as e:
        emit('script_error', {'error': str(e)}, namespace='/vcf')

@app.route('/bam_merger', methods=['GET', 'POST'])
@role_requis('superadmin') 
def bam_merger():
    if request.method == 'POST':
        input_dir = request.form['input_dir']
        output_dir = request.form['output_dir']
        
        if not all([input_dir, output_dir]):
            # Gestion d'erreur simplifiée pour l'exemple
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

@app.route('/generate_bam_script', methods=['GET'])
@role_requis('superadmin') 
def generate_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"samtools merge \"{config['output_dir']}/merged.bam\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"echo \"Merging complete for BAM files in {config['input_dir']}\"\n\n"
    return jsonify(script=script_content)

@app.route('/download_bam_script', methods=['GET'])
@role_requis('superadmin')
def download_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"samtools merge \"{config['output_dir']}/merged.bam\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"echo \"Merging complete for BAM files in {config['input_dir']}\"\n\n"
    
    # Utiliser un chemin absolu temporaire approprié
    script_path = '/data/Script_Site/tmp/bam_merge_script.sh'  # Assurez-vous que ce dossier existe et est accessible
    # script_path = 'C:/Users/aleks/OneDrive/Bureau/CHU-WebApp/tmp/bam_merge_script.sh'
    with open(script_path, 'w') as file:
        file.write(script_content)
    
    # Création de la réponse
    response = make_response(send_file(script_path, as_attachment=True, download_name="bam_merge_script.sh"))
    response.headers["Content-Disposition"] = "attachment; filename=bam_merge_script.sh"
    return response

    
@socketio.on('start_script', namespace='/test')
def handle_script():
    print("Received start script event from client.")
    emit("yoooo", namespace='/test')
    new_workflow = Workflow(name="BAM Merge", status="Running")
    db.session.add(new_workflow)
    db.session.commit()
    script_command = "bash /data/Script_Site/tmp/bam_merge_script.sh"
    try:
        process = subprocess.Popen(shlex.split(script_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        if stdout:
            for line in stdout.splitlines():
                emit('script_output', {'message': line, 'type': 'stdout'}, namespace='/test')
                new_workflow.status = "Completed"
                emit('script_error', {'status': "Completed"}, namespace='/test')

        if stderr:
            emit('script_output', {'message': stderr, 'type': 'stderr'}, namespace='/test')
            new_workflow.status = "Failed"
            emit('script_error', {'status': "Failed"}, namespace='/test')

        db.session.commit()

    except Exception as e:
        emit('script_error', {'error': str(e)}, namespace='/test')




@app.route('/basecalling', methods=['GET', 'POST'])
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
            return redirect(url_for('basecalling'))
        
        # Ajouter la configuration
        configurations_basecalling.append({
            "base_output_dir": base_output_dir,
            "input_dir": input_dir,
            "ref_genome": ref_genome_path,
            "qs_scores": qs_scores,
            "cuda_device": cuda_device,
            "model": model,
            "kit_name": kit_name
        })
        flash('Configuration added successfully.', 'success')
        return redirect(url_for('basecalling'))

    return render_template('index.html', configurations=configurations_basecalling)

@app.route('/generate_script', methods=['GET'])
@role_requis('superadmin') 
def generate_script():
    script_content = "#!/bin/bash\n\nsource ~/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_basecalling:
        qs_scores_list = config['qs_scores'].split()
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
${{DORADO_BIN}} basecaller -x "{config['cuda_device']}" --min-qscore "{qscore}" --no-trim --emit-fastq ${{MODEL_PATH}} ${{INPUT_DIR}} | \\
${{DORADO_BIN}} demux --kit-name "{config['kit_name']}" --emit-fastq --output-dir "${{OUTPUT_DIR}}"
echo "Processing complete for {config['input_dir']} with Q-score {qscore}"
"""
            script_content += f"for fastq_file in \"${{OUTPUT_DIR}}\"/*.fastq; do\n"
            script_content += f"    bam_file=\"${{fastq_file%.fastq}}.bam\"\n"
            script_content += f"    echo \"Aligning ${{fastq_file}} to reference genome...\"\n"
            script_content += f"    minimap2 -ax map-ont \"{config['ref_genome']}\" \"$fastq_file\" | samtools sort -o \"$bam_file\"\n"
            script_content += f"    samtools index \"$bam_file\"\n"
            script_content += f"    echo \"Alignment and BAM conversion completed for ${{bam_file}}\"\n"
            script_content += "done\n"
    script_content += "echo \"All processes are complete.\"\n"
    return jsonify(script=script_content)

@app.route('/add_config', methods=['POST'])
@role_requis('superadmin') 
def add_configuration():
    try:
        base_output_dir = request.form['base_output_dir']
        input_dir = request.form['input_dir']
        ref_genome_path = request.form['ref_genome']
        qs_scores = request.form['qs_scores']
        cuda_device = request.form['cuda_device']
        model = request.form['model']
        kit_name = request.form['kit_name']

        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.makedirs(app.config['UPLOAD_FOLDER'])
    

        new_config = {
            "base_output_dir": base_output_dir,
            "input_dir": input_dir,
            "ref_genome": ref_genome_path,
            "qs_scores": qs_scores,
            "cuda_device": cuda_device,
            "model": model,
            "kit_name": kit_name
        }
        configurations_basecalling.append(new_config)
        
        new_config_db = ConfigurationBasecalling(
            base_output_dir=base_output_dir,
            input_dir=input_dir,
            ref_genome=ref_genome_path,
            qs_scores=qs_scores,
            cuda_device=cuda_device,
            model=model,
            kit_name=kit_name
        )
        db.session.add(new_config_db)
        db.session.commit()
        
        return jsonify(success=True, configurations=configurations_basecalling)
    except Exception as e:
        return jsonify(success=False, message=str(e))
    
@app.route('/get_configurations_basecalling', methods=['GET'])
@role_requis('superadmin') 
def get_configurations_basecalling():
    return jsonify(configurations_basecalling)

@app.route('/get_configurations_merge', methods=['GET'])
@role_requis('superadmin') 
def get_configurations_merge():
    return jsonify(configurations_merge)

@app.route('/get_configurations_vcf', methods=['GET'])
@role_requis('superadmin') 
def get_configurations_vcf():
    return jsonify(configurations_vcf)

@app.route('/get_configurations_full_workflow', methods=['GET'])
@role_requis('superadmin') 
def get_configurations_full_workflow():
    return jsonify(configurations_full_workflow)


@app.route('/delete_config', methods=['POST'])
@role_requis('superadmin') 
def delete_configuration():
    index = request.json['index']
    try:
        configurations_basecalling.pop(index)
        return jsonify(success=True, configurations=configurations_basecalling)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/delete_config_merge', methods=['POST'])
@role_requis('superadmin') 
def delete_configuration_bam():
    index = request.json['index']
    try:
        configurations_merge.pop(index)
        return jsonify(success=True, configurations=configurations_merge)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/delete_config_vcf', methods=['POST'])
@role_requis('superadmin') 
def delete_configuration_vcf():
    index = request.json['index']
    try:
        configurations_vcf.pop(index)
        return jsonify(success=True, configurations=configurations_vcf)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/delete_config_full_workflow', methods=['POST'])
@role_requis('superadmin') 
def delete_configuration_full_workflow():
    index = request.json['index']
    try:
        configurations_full_workflow.pop(index)
        return jsonify(success=True, configurations=configurations_full_workflow)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/history-basecalling')
@role_requis('superadmin') 
def history():
    configurations = ConfigurationBasecalling.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-basecalling.html', configurations=configurations)

@app.route('/history-merge')
@role_requis('superadmin') 
def history_merge():
    configurations = ConfigurationMerge.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-merge.html', configurations=configurations)

@app.route('/history-vcf')
@role_requis('superadmin') 
def history_vcf():
    configurations = ConfigurationVCF.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-vcf.html', configurations=configurations)

@app.route('/history-vcf')
@role_requis('superadmin') 
def history_full_workflow():
    configurations = FullWorkflowConfiguration.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-full_workflow.html', configurations=configurations)

if __name__ == '__main__':
    with app.app_context():
        db.create_all()
    socketio.run(app, host='0.0.0.0',debug=True)
