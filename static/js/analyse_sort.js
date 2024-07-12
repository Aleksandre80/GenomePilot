function sortTable(column) {
    var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
    table = document.querySelector(".table");
    switching = true;
    // Définir la direction du tri comme ascendant:
    dir = "asc"; 
    // Faire une boucle qui continuera jusqu'à ce qu'aucun échange n'ait été fait:
    while (switching) {
        // Commencer par dire qu'aucun échange n'est fait:
        switching = false;
        rows = table.rows;
        // Boucler sur toutes les rangées du tableau (sauf la première, qui contient les en-têtes):
        for (i = 1; i < (rows.length - 1); i++) {
            // Commencer par dire qu'il ne devrait pas y avoir d'échange:
            shouldSwitch = false;
            // Obtenir les deux éléments à comparer, un de la ligne actuelle et un de la suivante:
            x = rows[i].getElementsByTagName("TD")[column];
            y = rows[i + 1].getElementsByTagName("TD")[column];
            // Vérifier si les deux rangées devraient changer de place, selon la direction, ascendant ou descendant:
            if (dir == "asc") {
                if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                    // Si oui, marquer comme un échange et casser la boucle:
                    shouldSwitch = true;
                    break;
                }
            } else if (dir == "desc") {
                if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
                    shouldSwitch = true;
                    break;
                }
            }
        }
        if (shouldSwitch) {
            // Si un échange a été marqué, faire l'échange et marquer que l'échange a été fait:
            rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
            switching = true;
            // Chaque fois qu'un échange est complété, augmenter ce compte:
            switchcount ++;      
        } else {
            // Si aucun échange n'a été fait et que la direction est "asc",
            // définir la direction à "desc" et refaire la boucle marquée:
            if (switchcount == 0 && dir == "asc") {
                dir = "desc";
                switching = true;
            }
        }
    }
}

document.addEventListener('DOMContentLoaded', function () {
    const checkboxes = document.querySelectorAll('.form-check-input');
    checkboxes.forEach(checkbox => {
        checkbox.addEventListener('change', filterTable);
    });

    function filterTable() {
        const activeFilters = Array.from(checkboxes)
            .filter(checkbox => checkbox.checked)
            .map(checkbox => checkbox.value.toLowerCase());

        const rows = document.querySelectorAll(".table tbody tr");
        rows.forEach(row => {
            const statusCell = row.querySelector("td:last-child span").textContent.toLowerCase();
            row.style.display = activeFilters.includes(statusCell) ? "" : "none";
        });
    }
});

document.addEventListener('DOMContentLoaded', function() {
    const searchInput = document.getElementById('searchInput');

    searchInput.addEventListener('keyup', function() {
        const searchValue = searchInput.value.toLowerCase();
        const rows = document.querySelectorAll(".table tbody tr");

        rows.forEach(row => {
            const nameCellText = row.cells[0].textContent.toLowerCase();
            if (nameCellText.includes(searchValue)) {
                row.style.display = "";
            } else {
                row.style.display = "none";
            }
        });
    });
});
