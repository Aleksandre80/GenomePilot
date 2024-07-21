from flask import session, flash, redirect, url_for
from functools import wraps

def role_requis(*roles_requis):
    def wrapper(fn):
        @wraps(fn)
        def decorated_view(*args, **kwargs):
            if get_role_utilisateur() not in roles_requis:
                flash("Accès refusé. Veuillez vous connecter avec un compte autorisé.", "warning")
                return redirect(url_for('common_bp.login'))
            return fn(*args, **kwargs)
        return decorated_view
    return wrapper

def get_role_utilisateur():
    return session.get('role', 'user')
