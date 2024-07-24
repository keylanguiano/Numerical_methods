let app = document.getElementById('typewriter');
 
let typewriter = new Typewriter(app, {
  loop: true,
  delay: 75,
  cursor: "<span style='color: #2dc4c6'>|</span>"
});
 
typewriter
  .pauseFor(2500)
  .typeString('<span style="color: #2dc4c6">Estudiante de Ingeniería en Sistemas Computacionales.</span>')
  .pauseFor(200)
  .deleteChars(10)
  .start();
