import os

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'votre_clé_secrète_ici'
    UPLOAD_FOLDER = 'path_to_save'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///configurations.db'
    SQLALCHEMY_TRACK_MODIFICATIONS = False
