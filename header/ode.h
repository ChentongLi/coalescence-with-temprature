typedef struct ODE
{
    double I;
    double lamb;
    double t;
}ODE;

double Tempreture(double t)
{   
    while(1){
        if (t<1) t+=365;
        else break;
    }
    int x=(int)(t-1);
    return TEM[x];
}

double birds(double t)
{
    if (t<0) t+=365;
    if (t<273 && t>74) return 0.0;
    if (t>0 && t<=74) return cos(M_PI/83*(t+92)+M_PI)+1;
    return cos(M_PI/83*(t-273)+M_PI)+1;
}

ODE *solution;
double dt=0.01;

double A=0.01;
double B=0.01;
double C=10;
double mu=10;
double Te=5;
double I0=5000;

ODE solver(ODE tmp){

    ODE result;
    result.t=tmp.t+dt;
    double tp=(Tempreture(tmp.t)-Te)/mu;
    result.I=tmp.I+dt*((A+B*exp(-tp*tp))*tmp.I+C*birds(tmp.t));
    result.lamb=1.0/result.I;
    return result;
}

void SolveOde(double T1_time, double End_time){

    int i;
    if(solution!=NULL) {
        free(solution);
        solution=NULL;
    }
    int n=ceil((End_time-T1_time)/dt)+1;
    solution=malloc(sizeof(ODE)*n);
    solution[0].I=I0;
    solution[0].lamb=1.0/I0;
    solution[0].t=T1_time;
    for(i=1;i<n;i++){
        solution[i]=solver(solution[i-1]);
    } 

  //  for(i=0;i<n;i++)
  //      printf("%e %e %e %e\n",solution[i].I,solution[i].t,solution[i].lamb);
}
