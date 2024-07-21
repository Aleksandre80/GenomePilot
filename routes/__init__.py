from flask import Blueprint

# Import des Blueprints individuels
from .basecalling import basecalling_bp
from .merge import merge_bp
from .vcf import vcf_bp
from .full_workflow import full_workflow_bp
from .anomalie_structure import anomalie_structure_bp
from .common import common_bp

# Fonction pour enregistrer les Blueprints dans l'application Flask
def register_blueprints(app):
    app.register_blueprint(basecalling_bp)
    app.register_blueprint(merge_bp)
    app.register_blueprint(vcf_bp)
    app.register_blueprint(full_workflow_bp)
    app.register_blueprint(anomalie_structure_bp)
    app.register_blueprint(common_bp)
