// VolRender.cpp : Defines the entry point for the console application.
//
#include <cstring>
#include <GL/glew.h>
#include <GL/glut.h>

#include "Raw3DData.h"
#include "transform.h"

#define PI 3.1415926535897932384626433832795

//#define PREINT

// GL ERROR CHECK
int CheckGLError(char *file, int line)
{
    GLenum glErr;
    int    retCode = 0;

    glErr = glGetError();
    while (glErr != GL_NO_ERROR) {
        printf("GL Error #%d ( %s ) in File %s at line: %d\n", glErr, gluErrorString(glErr),file, line );
        retCode = 1;
        glErr = glGetError();
    }
    return retCode;
}
#define CHECK_GL_ERROR() CheckGLError(__FILE__, __LINE__)

// DATA for EXAMPLE RENDERING

int W=1000;
int H=800;

Raw3DData obj;
unsigned char *tf;

float mm[16];
mat3<float> mat;
vec3<float> vi(1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0),
vj,
vk(-1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0));

int mov=0;
int posx, posy;
float alpha=0.0, beta=0.0;

unsigned int vshader;
const int HV_VOLDATA_NVS = 21;
const char *hvVolDataVS[HV_VOLDATA_NVS] = {
    "#version 140\n",
    "\n",
    "in vec4 vpos;\n",
    "uniform mat4 trans;\n",
    "uniform mat4 rev;\n",
    "uniform mat4 persp;\n",
    "uniform vec3 eye;\n",
    "uniform vec3 light;\n",
    "\n",
    "smooth out vec3 textcoord;\n",
    "noperspective out vec3 dir;\n",
    "noperspective out vec3 lidir;\n",
    //"smooth out vec3 view;\n",
    "\n",
    "void main(void) {\n",
    " vec4 v = trans * vpos;\n",
    " gl_Position = persp * v;\n",
    " vec3 view = v.xyz/v.w-eye;\n",
    " dir = normalize(mat3(rev) * view);\n",
    " lidir = normalize(mat3(rev) * light);\n",
    " textcoord = vec3((vpos.x+1.0)/2.0,(vpos.y+1.0)/2.0,(vpos.z+1.0)/2.0);\n",
    "}\n"
};

unsigned int fshader;
// USING NO PRE-INTEGRATION
const int HV_VOLDATA_NFS = 31;
const char *hvVolDataFS[HV_VOLDATA_NFS] = {
    "#version 140\n",
    "\n",
    "smooth in vec3 textcoord;\n",
    "noperspective in vec3 dir;\n",
    "out vec3 gl_FragColor;\n",
    "uniform sampler1D tf;\n",
    "uniform sampler3D data;\n",
    "noperspective in vec3 lidir;\n",
    "uniform vec4 lipara;\n",
    "uniform float slice, sampling;\n",
    "void main(void) {\n",
    " vec3 ndir = normalize(dir);\n",
    " vec3 ldir = normalize(lidir);\n",
    " vec3 dstep = ndir*(sampling/128.0);\n",
    " vec3 pos = textcoord+dstep;\n",
    " vec4 col=vec4(0.0,0.0,0.0,1.0);\n",
    " while(pos.x<1.0 && pos.x>0.0 && pos.y<1.0 && pos.y>0.0 && pos.z<1.0 && pos.z>0.0  && col.w>0.01)\n",
    "  {\n",
    "   vec4 val = vec4(texture(data, pos));\n",
    "   vec4 cc = vec4(texture(tf, val.w));\n",
    "   vec3 grad = normalize(vec3(val.x*2.0-1.0, val.y*2.0-1.0, val.z*2.0-1.0));\n",
    "	float alpha = col.w*cc.w*sampling;\n",
    "	vec3 half; if (dot(grad,ldir)<0.0 || dot(grad,ndir)>=0.0) half=vec3(0.0);\n",
    "	else half = normalize(ldir-ndir);\n",
    "	float dd = clamp(dot(ldir,grad),0.0,1.0); \n",
    "	float ds = pow(clamp(dot(grad,half),0.0,1.0),lipara.w);\n",
    "	col = vec4(col.xyz+(cc.xyz*(lipara.x+lipara.y*dd)+vec3(ds*lipara.z))*alpha,col.w*(1.0-cc.w*sampling));\n",
    "   pos = pos+dstep;\n",
    "  }\n",
    " gl_FragColor = col.xyz+vec3(col.w);\n",
    "}"
};

// USING CLASSIC  PRE-INTEGRATION
const int HV_VOLDATA_NFS2 = 35;
const char *hvVolDataFS2[HV_VOLDATA_NFS2] = {
    "#version 140\n",
    "\n",
    "smooth in vec3 textcoord;\n",
    "noperspective in vec3 dir;\n",
    "out vec3 gl_FragColor;\n",
    "uniform sampler2D pretab;\n",
    "uniform sampler3D data;\n",
    "noperspective in vec3 lidir;\n",
    "uniform vec4 lipara;\n",
    "\n",
    "void main(void) {\n",
    " vec3 ndir = normalize(dir);\n",
    " vec3 ldir = normalize(lidir); \n",
    " vec3 step = ndir*(0.5/128.0);\n",
    " vec3 pos = textcoord+step;\n",
    " vec4 col=vec4(0.0,0.0,0.0,1.0);\n",
    " vec4 oldval = vec4(texture(data, pos));\n",
    " vec3 oldgrad = normalize(vec3(oldval.x*2.0-1.0, oldval.y*2.0-1.0, oldval.z*2.0-1.0));\n",
    " pos = pos+step*0.5;\n",
    " while(pos.x<1.0 && pos.x>0.0 && pos.y<1.0 && pos.y>0.0 && pos.z<1.0 && pos.z>0.0 && col.w>0.01)\n",
    "  {\n",
    "   vec4 val = vec4(texture(data, pos));\n",
    "   vec4 cc = vec4(texture(pretab, vec2(oldval.w, val.w)));\n",
    "	float alpha = col.w;\n",
    "   vec3 grad = normalize(vec3(val.x*2.0-1.0, val.y*2.0-1.0, val.z*2.0-1.0));\n",
    "	vec3 half; if (dot(grad,ldir)<0.0 || dot(grad,ndir)>=0.0) half=vec3(0.0);\n",
    "	else half = normalize(ldir-ndir);\n",
    "	float dd = clamp(dot(ldir,grad),0.0,1.0); \n",
    "	float ds = clamp((dot(grad,half)-1.0)*lipara.w+1.0,0.0,1.0);\n",
    "	col = vec4(col.xyz+(cc.xyz*(lipara.x+lipara.y*dd)+vec3(ds*lipara.z))*alpha,alpha*(1.0-cc.w));\n",
    "   pos = pos+step;\n",
    "	oldval = val; oldgrad=grad; \n",
    "  }\n",
    " gl_FragColor = col.xyz+vec3(col.w);\n",
    "}"
};

unsigned int id[4], textid[4];
unsigned int gpuprog;
int itrans, iposition, ipersp, ilight, ieye, ipara, islice, isampling;
int isampler, itf;

float vert[3*8] = { -1.0f, -1.0f, -1.0f,   -1.0f, -1.0f, 1.0f,   -1.0f, 1.0f, -1.0f,   -1.0f, 1.0f, 1.0f,
                    1.0f, -1.0f, -1.0f,   1.0f, -1.0f, 1.0f,   1.0f, 1.0f, -1.0f,   1.0f, 1.0f, 1.0f };
//float vert[3*4] = { 0.0f,-1.0f,-1.0f,0.0f,-1.0f,1.0f, 0.0f,1.0f,1.0f, 0.0f,1.0f,-1.0f };
unsigned int faces[36] = { 0,1,2, 2,1,3 , 3,1,5, 3,5,7,   0,2,4, 4,2,6, 1,0,4, 1,4,5, 2,3,7, 2,7,6, 6,7,5, 5,4,6 };

void createlist()
{
    int i;
    int compiled=0;
    char buffer[1024];

    // create buffers for vertices, normals and element indices
    glGenBuffers(4,id);
    glBindBuffer(GL_ARRAY_BUFFER, id[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*8*sizeof (float), vert, GL_STATIC_DRAW);
    //CHECK_GL_ERROR();

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*12*sizeof (unsigned int), faces, GL_STATIC_DRAW);
    //CHECK_GL_ERROR();

    // vertex shader
    int vshaderlength[500];
    vshader = glCreateShader(GL_VERTEX_SHADER);
    for (i=0; i<HV_VOLDATA_NVS; i++) vshaderlength[i]=strlen(hvVolDataVS[i]);
    glShaderSource(vshader, HV_VOLDATA_NVS, (const GLchar **)hvVolDataVS, vshaderlength);
    glCompileShader(vshader);
    //CHECK_GL_ERROR();
    glGetShaderiv(vshader, GL_COMPILE_STATUS, &compiled);
    //CHECK_GL_ERROR();
    if (compiled) { printf("Vertex Shader compilation successful!\n"); }
    else
    {
        printf("Shader compile error:\n");
        glGetShaderInfoLog(vshader, 1024, &i, buffer);
        printf("%s\n", buffer);
    }
    // fragment shader
    int fshaderlength[500];
    fshader = glCreateShader(GL_FRAGMENT_SHADER);
#ifdef PREINT
    for (i=0; i<HV_VOLDATA_NFS2; i++) fshaderlength[i]=strlen(hvVolDataFS2[i]);
    glShaderSource(fshader, HV_VOLDATA_NFS2, (const GLchar **)hvVolDataFS2, fshaderlength);
#else
    for (i=0; i<HV_VOLDATA_NFS; i++) fshaderlength[i]=strlen(hvVolDataFS[i]);
    glShaderSource(fshader, HV_VOLDATA_NFS, (const GLchar **)hvVolDataFS, fshaderlength);
#endif
    glCompileShader(fshader);
    //CHECK_GL_ERROR();
    glGetShaderiv(fshader, GL_COMPILE_STATUS, &compiled);
    //CHECK_GL_ERROR();
    if (compiled) { printf("Fragment Shader compilation successful!\n"); }
    else
    {
        printf("Shader compile error:\n");
        glGetShaderInfoLog(fshader, 1024, &i, buffer);
        printf("%s\n", buffer);
    }

    gpuprog = glCreateProgram();
    glAttachShader(gpuprog, vshader);
    glAttachShader(gpuprog, fshader);
    glLinkProgram(gpuprog);
    //CHECK_GL_ERROR();
    glGetProgramiv(gpuprog, GL_LINK_STATUS, &compiled);
    //CHECK_GL_ERROR();
    if (compiled) { printf("Link successful!\n"); }
    else
    {
        printf("Link error:\n");
        glGetProgramInfoLog(gpuprog, 1024, &i, buffer);
        printf("%s\n", buffer);
    }

    glUseProgram(gpuprog);
    //CHECK_GL_ERROR();

    int nattr; glGetProgramiv(gpuprog, GL_ACTIVE_ATTRIBUTES, &nattr);
    printf("activated attributes : %d\n", nattr);
    int nunif; glGetProgramiv(gpuprog, GL_ACTIVE_UNIFORMS, &nunif);
    printf("activated uniforms : %d\n", nunif);

    int ivpos = glGetAttribLocation(gpuprog, "vpos");
    printf("attribute vpos = %d\n", ivpos);
    //CHECK_GL_ERROR();

    itrans = glGetUniformLocation(gpuprog, "trans");
    //CHECK_GL_ERROR();
    printf("uniform trans = %d\n", itrans);
    iposition = glGetUniformLocation(gpuprog, "rev");
    //CHECK_GL_ERROR();
    printf("uniform rev = %d\n", iposition);
    ipersp = glGetUniformLocation(gpuprog, "persp");
    //CHECK_GL_ERROR();
    printf("uniform persp = %d\n", ipersp);
    ieye = glGetUniformLocation(gpuprog, "eye");
    //CHECK_GL_ERROR();
    printf("uniform eyepos = %d\n", ieye);
    ilight = glGetUniformLocation(gpuprog, "light");
    //CHECK_GL_ERROR();
    printf("uniform light = %d\n", ilight);
    ipara = glGetUniformLocation(gpuprog, "lipara");
    //CHECK_GL_ERROR();
    printf("uniform light = %d\n", ilight);
    isampler = glGetUniformLocation(gpuprog, "data");
    //CHECK_GL_ERROR();
    printf("sampler3D data = %d\n", isampler);

#ifdef PREINT
    itf = glGetUniformLocation(gpuprog, "pretab");
    //CHECK_GL_ERROR();
    printf("sampler1D pretab = %d\n", itf);
#else
    islice = glGetUniformLocation(gpuprog, "slice");
    //CHECK_GL_ERROR();
    printf("uniform slice = %d\n", islice);
    isampling = glGetUniformLocation(gpuprog, "sampling");
    //CHECK_GL_ERROR();
    printf("uniform sampling = %d\n", isampling);
    itf = glGetUniformLocation(gpuprog, "tf");
    //CHECK_GL_ERROR();
    printf("sampler1D tf = %d\n", itf);
#endif

    glBindBuffer(GL_ARRAY_BUFFER, id[0]);
    glVertexAttribPointer(ivpos, 3, GL_FLOAT, false, 0, NULL);
    //CHECK_GL_ERROR();
    glEnableVertexAttribArray(ivpos);
    //CHECK_GL_ERROR();
}

void inittextures()
{ 
    glGenTextures(4, textid);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textid[0]);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, obj.getSx(), obj.getSy(), obj.getSz(), 0, GL_RGBA, GL_UNSIGNED_BYTE, obj.getGrad());
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

#ifdef PREINT
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_2D, textid[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, obj.getTabSize(), obj.getTabSize(), 0, GL_RGBA, GL_UNSIGNED_SHORT, obj.getPreIntTab());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
#else
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_1D, textid[1]);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA8, obj.getTabSize(), 0, GL_RGBA, GL_UNSIGNED_BYTE, obj.getTransfer());
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
#endif
}

// OPENGL DISPLAY
void display()
{
    transform<float> tr;

    // effacer l'\E9cran
    glClearColor(1.0,1.0,1.0,0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

    // initialiser le viewport
    glViewport(0, 0, W, H);

    // activer le Z-buffer
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glFrontFace(GL_CCW);
    glEnable(GL_CULL_FACE);


    float li[3] = { 1.0, 1.0, 1.0 };
    glUniform3fv(ilight, 1, li);
    //CHECK_GL_ERROR();
    float liparadata[4];
    // define Ka, Kd,Ks, Phong exp
    liparadata[0]=0.3; liparadata[1]=0.7; liparadata[2]=0.1; liparadata[3]=60.0;
    glUniform4fv(ipara, 1, liparadata);
    //CHECK_GL_ERROR();

#ifdef PREINT
#else
    float slice_val=0.0;
    glUniform1fv(islice, 1, &slice_val);
    float sampling_val= 0.5 ; //voxelstep;
    glUniform1fv(isampling, 1, &sampling_val);
#endif
    // initialiser la transformation perspective
    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glFrustum(-0.42, 0.42, -0.42*H/W, 0.42*H/W, 1.0, 50.0);


    // initialiser la transformation de l'objet
    mat = mat3<float>(vi,vj,vk);
    mat.inverse();
    mm[0]=mat.I().X(); mm[1]=mat.I().Y(); mm[2]=mat.I().Z(); mm[3]=0.0;
    mm[4]=mat.J().X(); mm[5]=mat.J().Y(); mm[6]=mat.J().Z(); mm[7]=0.0;
    mm[8]=mat.K().X(); mm[9]=mat.K().Y(); mm[10]=mat.K().Z(); mm[11]=0.0;
    mm[12]=0.0; mm[13]=0.0; mm[14]=0.0; mm[15]=1.0;

    float eye[3]={ -12.0, -12.0, 12.0 };
    glUniform3fv(ieye, 1, eye);
    //CHECK_GL_ERROR();

    tr.setIdentity();
    tr.frustum(-0.1, 0.1, -0.1*H/W, 0.1*H/W, 0.08, 50.0);
    tr.frame(vi,vj,vk);
    tr.translation(vec3<float>(-eye[0], -eye[1], -eye[2]));
    tr.getMatrix(mm);
    glUniformMatrix4fv(ipersp, 1, false, mm);
    //CHECK_GL_ERROR();

    tr.setIdentity();
    tr.rotation(vec3<float>(0.0,0.0,1.0), alpha/180.0*PI);
    tr.rotation(vec3<float>(1.0,0.0,0.0), beta/180.0*PI);
    tr.scale(vec3<float>(8.0,8.0,8.0));
    tr.getMatrix(mm);
    glUniformMatrix4fv(itrans, 1, false, mm);
    //CHECK_GL_ERROR();
    tr.setIdentity();
    tr.scale(vec3<float>(1.0/8.0,1.0/8.0,1.0/8.0));
    tr.rotation(vec3<float>(1.0,0.0,0.0), -beta/180.0*PI);
    tr.rotation(vec3<float>(0.0,0.0,1.0), -alpha/180.0*PI);
    tr.getMatrix(mm);
    glUniformMatrix4fv(iposition, 1, false, mm);
    //CHECK_GL_ERROR();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textid[0]);
    //CHECK_GL_ERROR();
    glUniform1i(isampler,0);
    //CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_3D);
    //CHECK_GL_ERROR();
    glActiveTexture(GL_TEXTURE0+1);
#ifdef PREINT
    glBindTexture(GL_TEXTURE_2D, textid[1]);
    //CHECK_GL_ERROR();
    glUniform1i(itf,1);
    //CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_2D);
    //CHECK_GL_ERROR();
#else
    glBindTexture(GL_TEXTURE_1D, textid[1]);
    //CHECK_GL_ERROR();
    glUniform1i(itf,1);
    //CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_1D);
    //CHECK_GL_ERROR();
#endif

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id[1]);
    glDrawElements(GL_TRIANGLES, 12*3, GL_UNSIGNED_INT, NULL);
    //CHECK_GL_ERROR();

    glutSwapBuffers();
}


// MANAGE EVENTS
void clic(int button, int state, int x, int y)
{
    if (button==0 && state==0) { mov=1; posx=x; posy=y; }
    else mov=0;
}

void bouge(int x, int y)
{
    if (mov==1)
    {
        //printf("move %d,%d\n", x,y);
        alpha += (float)(x-posx)*0.5;
        posx=x;
        beta += (float)(y-posy)*0.25;
        posy=y;
        display();
    }
}

void key(unsigned char c, int x, int y)
{
    printf("char=%c (%d)\n",c, (int)c);
}



// START PROGRAMM AND INITIALIZE WINDOW
int main(int argc, char **argv)
{
    vj.cross(vk,vi);
    //vj.reverse();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowPosition(0,0) ;
    glutInitWindowSize(W, H);
    glutCreateWindow("ma fenetre");

    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }

    glutDisplayFunc(display);
    glutMouseFunc(clic);
    glutMotionFunc(bouge);
    //glutKeyboardFunc(key);

    // fonction de transfert
    tf = new unsigned char [128*4];
    int i;
    for (i=0; i<127; i++) { tf[4*i]=0; tf[4*i+1]=0; tf[4*i+2]=0; tf[4*i+3]=0; }
    for (i=10; i<30; i++) { tf[4*i]=static_cast<unsigned char>(0.99*255.0); tf[4*i+1]=(unsigned char)(0.9*255.0); tf[4*i+2]=(unsigned char)(0.6*255.0); tf[4*i+3]=(unsigned char)(0.025*255.0); }
    for (i=30; i<40; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.29*255.0); tf[4*i+2]=(unsigned char)(0.1*255.0); tf[4*i+3]=(unsigned char)(0.2*255.0); }
    for (i=40; i<50; i++) { tf[4*i]=(unsigned char)(0.9*255.0); tf[4*i+1]=(unsigned char)(0.7*255.0); tf[4*i+2]=(unsigned char)(0.9*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }
    for (i=70; i<90; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.8*255.0); tf[4*i+2]=(unsigned char)(0.7*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }
    for (i=100; i<126; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.9*255.0); tf[4*i+2]=(unsigned char)(0.6*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }

    // charger les donnees
    printf("loading raw data...\n");
    obj.importRAW("data/CTKnee.raw",379,229,305, 1.0/255.0);
    //obj.computeGrad();
    obj.importGrad("data/CTKnee_grad.raw");
    obj.computePreIntegration(tf, 128, 0.5, 20);

    // init GPU
    createlist();
    inittextures();

    glutMainLoop();
    return 0;
}
