from flask import Blueprint

# Import des Blueprints individuels
from .basecalling import basecalling_bp
from .merge import merge_bp
from .vcf import vcf_bp
from .full_workflow import full_workflow_bp
from .anomalie_structure import anomalie_structure_bp
from .moyenne_reads import moyenne_reads_bp
from .phred_quality import phred_quality_bp
from .cible_reads import cible_reads_bp
from .reads_length import reads_length_bp
from .methylation import methylation_bp
from .methylartist import methylartist_bp
from .methylartist_violon import methylartist_violon_bp
from .common import common_bp

# Fonction pour enregistrer les Blueprints dans l'application Flask
def register_blueprints(app):
    app.register_blueprint(basecalling_bp)
    app.register_blueprint(merge_bp)
    app.register_blueprint(vcf_bp)
    app.register_blueprint(full_workflow_bp)
    app.register_blueprint(anomalie_structure_bp)
    app.register_blueprint(common_bp)
    app.register_blueprint(moyenne_reads_bp)
    app.register_blueprint(phred_quality_bp)
    app.register_blueprint(cible_reads_bp)
    app.register_blueprint(reads_length_bp)
    app.register_blueprint(methylation_bp)
    app.register_blueprint(methylartist_bp)
    app.register_blueprint(methylartist_violon_bp)
    
