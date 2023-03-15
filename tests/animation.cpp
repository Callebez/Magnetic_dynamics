#include <GL/glut.h>
#include <cmath>

// global variable for rotation angle
float angle = 0.0f;

void display()
{
   glClear(GL_COLOR_BUFFER_BIT); 
  glColor3f(1.0f, 1.0f, 1.0f);

  // draw circle
  glBegin(GL_TRIANGLE_FAN);
  
    // center of circle
    glVertex3f(0.0f, 0.0f, 0.0f);
    
    // half white
    glColor3f(1.0f, 1.0f, 1.0f);
    for (float i = 0.0f; i < M_PI+angle; i += 0.01f)
    {
        glVertex3f(cos(i), sin(i), 0.0f);
    }

  
  glEnd();

  glColor3f(0.5f, 0.5f, 0.5f);
  glBegin(GL_TRIANGLE_FAN);
    // half grey
    glVertex3f(0.0f, 0.0f, 0.0f);

    for (float i = M_PI + angle ; i <= 2.0f * M_PI + angle; i += 0.01f)
    {
        glVertex3f(cos(i), sin(i), 0.0f);
    }

  glEnd();

  // restore matrix state
  glPopMatrix();
  
  // increase angle for next frame
  // angle += 1.0f;
  
  // swap buffers to display the result
  glutSwapBuffers();
}
void updateCircle()
{
  angle+= 0.01f;
}
void idle() 
{
	// count=0;
  updateCircle();
	// std::cout<<"tempo: "<<fun.solver->rk_int->data[count].second<<"\n";
	// std::cout<< count<<"Valor: "<<fun.solver->rk_int->data[count].first[0]<<'\n';

    glutPostRedisplay();
}

int main(int argc, char** argv)
{
  // initialize GLUT
  // glutInit(&argc, argv);
  // glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  // glutInitWindowSize(640, 480);
  
  // // create window
  // glutCreateWindow("Spinning Half-White Half-Grey Circle");
  
  // // set display function
  // glutDisplayFunc(display);
  
  // // enable depth testing
  // glEnable(GL_DEPTH_TEST);
  
  // // enter main loop
  // glutMainLoop();

  	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(640, 480);
	glutCreateWindow("Three Circles with Arrows");
  glutDisplayFunc(display);
  // glutIdleFunc(idle);
    glutIdleFunc(idle);
  glutMainLoop();
  
  return 0;
}
