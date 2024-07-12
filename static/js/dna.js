let amplitude = 20; // Hauteur réduite de l'onde
let xspacing = 8; // Diminution de la distance entre chaque position horizontale
let dx; // Valeur pour l'incrément de x
let angle = 0; // Angle de départ
let w; // Largeur de l'onde entière
let period = 200; // Nombre de pixels avant que l'onde se répète
let yvalues, yvalues1; // Tableaux pour stocker les valeurs y

function setup(){
    createCanvas(150, 150); // Taille fixe plus petite pour le logo
    w = width * 0.8; // Adaptation de la largeur de l'onde à la largeur du canvas
    dx = (TWO_PI / period) * xspacing;
    yvalues = new Array(floor(w / xspacing));
    yvalues1 = new Array(floor(w / xspacing));
}

function draw(){
    clear(); // Utilisez clear() pour garder le fond transparent
    noFill();
    stroke('#aaa'); // Couleur du cadre
    strokeWeight(2); // Épaisseur du cadre
    ellipse(width / 2, height / 2, 140, 140); // Dessin du cercle cadre
    translate(width / 2, height / 2); // Centre le point de rotation
    rotate(PI / 4); // Inclinaison de 45 degrés
    translate(-width / 2, -height / 2); // Revenir au point d'origine
    calculate();
    update();
}

function calculate(){
    angle += 0.02;
    let x = angle;
    for(let i = 0; i < yvalues.length; i++){
        yvalues[i] = sin(x*2) * amplitude;
        yvalues1[i] = sin(2*PI - x*2) * amplitude;
        x += dx;
    }
}

function update(){
    translate(0, 40); // Ajustement de la position pour le petit canvas
    for(let x = 0; x < yvalues.length; x++){
        fill('#f70');
        stroke('#f70');
        strokeWeight(1.5); // Réduction de l'épaisseur du trait
        ellipse(width/2 + yvalues[x], x * xspacing, 8, 8); // Réduction de la taille des ellipses
        line(width/2, x * xspacing, width/2 + yvalues[x], x * xspacing);
    }
    for(let x = 0; x < yvalues1.length; x++){
        fill('#17f');
        stroke('#17f');
        strokeWeight(1.5);
        ellipse(width/2 + yvalues1[x], x * xspacing, 8, 8);
        line(width/2, x * xspacing, width/2 + yvalues1[x], x * xspacing);
    }
}

function windowResized(){
    resizeCanvas(150, 150); // Maintenir la taille fixe même lors du redimensionnement de la fenêtre
}
