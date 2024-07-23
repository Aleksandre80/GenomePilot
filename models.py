from datetime import datetime
from extensions import db

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
    output_dir = db.Column(db.String(120), nullable=False)
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

class ConfigurationAnomalieStructure(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_dir = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationMoyenneReads(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_file = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    bed_file = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationPhredQuality(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_file = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    chr = db.Column(db.String(120), nullable=False)
    pos1 = db.Column(db.String(120), nullable=False)
    pos2 = db.Column(db.String(120), nullable=False)
    phred_min = db.Column(db.String(120), nullable=False)
    logs = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationCibleReads(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_file = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    bed_file = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)
    
class ConfigurationReadsLength(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    input_file = db.Column(db.String(120), nullable=False)
    output_dir = db.Column(db.String(120), nullable=False)
    length_option = db.Column(db.String(120), nullable=False)
    min_length = db.Column(db.String(120), nullable=True)
    max_length = db.Column(db.String(120), nullable=True)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)

class Workflow(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False)
    launch_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    status = db.Column(db.String(20), nullable=False, default='Pending')

    def __repr__(self):
        return f"Workflow('{self.name}', '{self.launch_date}', '{self.status}')"