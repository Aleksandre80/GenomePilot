document.addEventListener('DOMContentLoaded', function () {
    const updateRunningCount = () => {
        fetch('/running_workflows_count')
            .then(response => response.json())
            .then(data => {
                document.getElementById('runningCount').textContent = data.running_workflows_count;
            })
            .catch(error => console.error('Erreur lors de la récupération des données:', error));
    };

    // Mise à jour périodique du nombre d'analyses en cours
    updateRunningCount();
    setInterval(updateRunningCount, 30000); // Met à jour toutes les 30 secondes
});