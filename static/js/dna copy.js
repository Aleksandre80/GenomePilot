let amplitude = 10; // Height of wave
let xspacing = 8; // Distance between each horizontal location
let dx; // Value for incrementing x
let angle = 0; // Start angle
let w; // Width of entire wave
let period = 200; // How many pixels before the wave repeats
let yvalues, yvalues1; // Array to store yvalues

function setup(){
    createCanvas(100, 100); // Réduire la taille du canvas pour un logo
    w = width; // Adapter la largeur de l'onde à la largeur du canvas
    dx = (TWO_PI / period) * xspacing;
    yvalues = new Array(floor(w / xspacing));
    yvalues1 = new Array(floor(w / xspacing));
}
function draw(){
    background(255); // Utiliser un fond clair pour le logo
    translate(width / 2, height / 2); // Centrer l'animation
    rotate(frameCount / 100.0); // Rotation continue pour l'effet de chargement
    translate(-width / 2, -height / 2);
    calculate();
    update();
}


function calculate(){
	angle += 0.02;
	let x = angle;
	for(let i = 0; i < yvalues.length; i += 1){
		yvalues[i] = sin(x*2) * amplitude;
		yvalues1[i] = sin(2*PI - x*2) * amplitude;
		x += dx;
	}
}
function update(){
    for(let x = 0; x < yvalues.length; x++){
        stroke('#f70');
        strokeWeight(2);
        line(width / 2, x * xspacing, width / 2 + yvalues[x], x * xspacing);
    }
    for(let x = 0; x < yvalues1.length; x++){
        stroke('#17f');
        strokeWeight(2);
        line(width / 2, x * xspacing, width / 2 + yvalues1[x], x * xspacing);
    }
}

function windowResized(){
  resizeCanvas(windowWidth, windowHeight);
}