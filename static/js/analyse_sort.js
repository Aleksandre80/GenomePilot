function filterTable() {
    const nameFilter = document.getElementById('nameFilter').value.toLowerCase();
    const table = document.getElementById('workflowTableBody');
    const rows = table.getElementsByTagName('tr');
    
    const pendingChecked = document.getElementById('filterPending').checked;
    const runningChecked = document.getElementById('filterRunning').checked;
    const completedChecked = document.getElementById('filterCompleted').checked;
    const failedChecked = document.getElementById('filterFailed').checked;
    
    for (let i = 0; i < rows.length; i++) {
        const cells = rows[i].getElementsByTagName('td');
        const nameCell = cells[0].textContent.toLowerCase().trim();
        const statusCell = cells[4].textContent.toLowerCase().trim();
        
        const nameMatches = nameCell.includes(nameFilter);
        const statusMatches = 
            ((statusCell === 'pending' && pendingChecked) ||
             (statusCell === 'running' && runningChecked) ||
             (statusCell === 'completed' && completedChecked) ||
             (statusCell === 'failed' && failedChecked));
        
        if (nameMatches && statusMatches) {
            rows[i].style.display = '';
        } else {
            rows[i].style.display = 'none';
        }
    }
}