from flask import Blueprint, render_template, request, jsonify, session, redirect, url_for, flash
from models import Workflow
from utils import role_requis, get_role_utilisateur
import psutil
from extensions import db
from datetime import datetime

common_bp = Blueprint('common_bp', __name__)

@common_bp.route('/')
@role_requis('superadmin') 
def accueil():
    return render_template('workflow.html')

@common_bp.route('/status', methods=['GET', 'POST'])
@role_requis('superadmin')
def status():
    workflows = Workflow.query.order_by(Workflow.start_time.desc()).all()
    # new_workflow = Workflow(name="BAM Merge", status="Running", start_time=datetime.utcnow(), output_dir="output_dir", end_time=datetime.utcnow())
    # db.session.add(new_workflow)
    # db.session.commit()
    return render_template('status.html', workflows=workflows)

@common_bp.route('/infos')
def infos():
    return render_template('infos.html')
    

@common_bp.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        password = request.form['password']
        username = request.form['username']
        role = "user"
        if password == "Aleksandre" and username == "admin":
            role = "superadmin"
        session['password'] = password
        session['role'] = role
        
        if role == "admin" or role == "superadmin":
            return redirect(url_for('common_bp.accueil'))
        else:
            return redirect(url_for('common_bp.login'))
            
    return render_template('login.html')

@common_bp.route('/diskinfo')
def disk_info():
    partitions = psutil.disk_partitions()
    disk_data = []
    for partition in partitions:
        try:
            usage = psutil.disk_usage(partition.mountpoint)
            disk_info = {
                'device': partition.device,
                'mountpoint': partition.mountpoint,
                'used': usage.used / (1024**3),
                'free': usage.free / (1024**3),
                'total': usage.total / (1024**3)
            }
            disk_data.append(disk_info)
        except OSError as e:
            disk_info = {
                'device': partition.device,
                'mountpoint': partition.mountpoint,
                'error': str(e)
            }
            disk_data.append(disk_info)
    return jsonify(disk_data)

@common_bp.route('/disk')
def disk():
    return render_template('disk.html')

@common_bp.route('/running_workflows_count')
def running_workflows_count():
    count = Workflow.query.filter_by(status='Running').count()
    return jsonify({'running_workflows_count': count})

@common_bp.route('/methylation-hub')
def methylation():
    return render_template('methylation-hub.html')

@common_bp.route('/hierarchy')
def hierarchy():
    return render_template('templates-pipelines.html')