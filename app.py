import os
from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

app = Flask(__name__)
app.secret_key = 'votre_clé_secrète_ici'
app.config['UPLOAD_FOLDER'] = 'path_to_save'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///configurations.db'
db = SQLAlchemy(app)

configurations_basecalling = []
configurations_merge = []
configurations_vcf = []

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


@app.route('/', methods=['GET', 'POST'])
def accueil():
    return render_template('home.html')
    
@app.route('/vcf_creator', methods=['GET', 'POST'])
def vcf_creator():
    if request.method == 'POST':
        # Extract file name and save the file
        ref_genome_file = request.files['ref_genome']
        bam_file = request.form['bam_file']
        output_vcf = request.form['output_vcf']

        if not all([ref_genome_file.filename, bam_file, output_vcf]):
            return jsonify(success=False, message="Please specify all fields.")

        # Create the save directory if it does not exist
        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.makedirs(app.config['UPLOAD_FOLDER'])

        ref_genome_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(ref_genome_file.filename))
        ref_genome_file.save(ref_genome_path)

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
def generate_vcf_script():
    script_content = "#!/bin/bash\n\nsource ~/miniconda3/etc/profile.d/conda.sh\nconda activate genomics\n\n"
    for config in configurations_vcf:
        script_content += f"samtools faidx {config['ref_genome']}\n"
        script_content += f"samtools index {config['bam_file']}\n"
        script_content += f"bcftools mpileup -Ou -f {config['ref_genome']} {config['bam_file']} | bcftools call -mv -Ob -o {config['output_vcf']}.bcf\n"
        script_content += "echo \"Variant calling and file processing completed.\"\n"
    return jsonify(script=script_content)


@app.route('/bam_merger', methods=['GET', 'POST'])
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
def generate_bam_script():
    script_content = "#!/bin/bash\n\n"
    for config in configurations_merge:
        script_content += f"mkdir -p \"{config['output_dir']}\"\n"
        script_content += f"samtools merge \"{config['output_dir']}/merged.bam\" \"{config['input_dir']}\"/*.bam\n"
        script_content += f"echo \"Merging complete for BAM files in {config['input_dir']}\"\n\n"
    return jsonify(script=script_content)


@app.route('/basecalling', methods=['GET', 'POST'])
def basecalling():
    if request.method == 'POST':
        base_output_dir = request.form['base_output_dir']
        input_dir = request.form['input_dir']
        ref_genome = request.files['ref_genome']
        qs_scores = request.form['qs_scores']
        cuda_device = request.form['cuda_device']
        model = request.form['model']
        kit_name = request.form['kit_name']

        if not all([base_output_dir, input_dir, ref_genome, qs_scores, cuda_device, model, kit_name]):
            flash('Please fill all fields before adding a configuration.', 'error')
            return redirect(url_for('basecalling'))
        
        # Créer le dossier de sauvegarde s'il n'existe pas déjà
        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.makedirs(app.config['UPLOAD_FOLDER'])
        
        ref_genome_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(ref_genome.filename))
        ref_genome.save(ref_genome_path)

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
${{DORADO_BIN}} basecaller -x "{config['cuda_device']}" --min-qscore "{qscore}" --no-trim --emit_fastq ${{MODEL_PATH}} ${{INPUT_DIR}} | \\
${{DORADO_BIN}} demux --kit-name "{config['kit_name']}" --emit_fastq --output-dir "${{OUTPUT_DIR}}"
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
def add_configuration():
    try:
        base_output_dir = request.form['base_output_dir']
        input_dir = request.form['input_dir']
        ref_genome_file = request.files['ref_genome']
        qs_scores = request.form['qs_scores']
        cuda_device = request.form['cuda_device']
        model = request.form['model']
        kit_name = request.form['kit_name']

        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.makedirs(app.config['UPLOAD_FOLDER'])
        
        ref_genome_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(ref_genome_file.filename))
        ref_genome_file.save(ref_genome_path)

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
def get_configurations_basecalling():
    return jsonify(configurations_basecalling)

@app.route('/get_configurations_merge', methods=['GET'])
def get_configurations_merge():
    return jsonify(configurations_merge)

@app.route('/get_configurations_vcf', methods=['GET'])
def get_configurations_vcf():
    return jsonify(configurations_vcf)


@app.route('/delete_config', methods=['POST'])
def delete_configuration():
    index = request.json['index']
    try:
        configurations_basecalling.pop(index)
        return jsonify(success=True, configurations=configurations_basecalling)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/delete_config_merge', methods=['POST'])
def delete_configuration_bam():
    index = request.json['index']
    try:
        configurations_merge.pop(index)
        return jsonify(success=True, configurations=configurations_merge)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/delete_config_vcf', methods=['POST'])
def delete_configuration_vcf():
    index = request.json['index']
    try:
        configurations_vcf.pop(index)
        return jsonify(success=True, configurations=configurations_vcf)
    except IndexError:
        return jsonify(success=False, message="Configuration not found")
    
@app.route('/history-basecalling')
def history():
    configurations = ConfigurationBasecalling.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-basecalling.html', configurations=configurations)

@app.route('/history-merge')
def history_merge():
    configurations = ConfigurationMerge.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-merge.html', configurations=configurations)

@app.route('/history-vcf')
def history_vcf():
    configurations = ConfigurationVCF.query.all()
    configurations.sort(key=lambda x: x.date_created, reverse=True)
    return render_template('history-vcf.html', configurations=configurations)

if __name__ == '__main__':
    with app.app_context():
        db.create_all()
    app.run(debug=True, host='0.0.0.0')
