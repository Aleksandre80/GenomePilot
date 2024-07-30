import os
from flask import Flask
from config import Config
from extensions import db, migrate, socketio

app = Flask(__name__)
app.config.from_object(Config)
db.init_app(app)
migrate.init_app(app, db)
socketio.init_app(app)

from routes.basecalling import basecalling_bp
from routes.merge import merge_bp
from routes.vcf import vcf_bp
from routes.full_workflow import full_workflow_bp
from routes.anomalie_structure import anomalie_structure_bp
from routes.moyenne_reads import moyenne_reads_bp
from routes.phred_quality import phred_quality_bp
from routes.cible_reads import cible_reads_bp
from routes.reads_length import reads_length_bp
from routes.methylation import methylation_bp
from routes.methylartist import methylartist_bp
from routes.methylartist_violon import methylartist_violon_bp
from routes.common import common_bp

app.register_blueprint(basecalling_bp)
app.register_blueprint(merge_bp)
app.register_blueprint(vcf_bp)
app.register_blueprint(full_workflow_bp)
app.register_blueprint(anomalie_structure_bp)
app.register_blueprint(moyenne_reads_bp)
app.register_blueprint(phred_quality_bp)
app.register_blueprint(cible_reads_bp)
app.register_blueprint(reads_length_bp)
app.register_blueprint(methylation_bp)
app.register_blueprint(methylartist_bp)
app.register_blueprint(methylartist_violon_bp)
app.register_blueprint(common_bp)


if __name__ == '__main__':
    with app.app_context():
        db.create_all()
    socketio.run(app, host='0.0.0.0', debug=True)
