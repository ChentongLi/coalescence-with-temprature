typedef struct ODE
{
    double I;
    double lamb;
    double t;
}ODE;

double Tempreture(double t,int city)
{   
    while(1){
        if (t<1) t+=365;
        else break;
    }
    int x=(int)(t-1);
    if(city==0) return TEM0[x];
    if(city==1) return TEM1[x];
}

ODE *solution[2];
double dt=0.01;

double A[2]={0.01,0.01};
double B[2]={0.01,0.01};
double C=10;
double mu[2]={10,10};
double Te=5;
double I0[2]={10000,10000};

ODE solver(ODE tmp,int i){

    ODE result;
    result.t=tmp.t+dt;
    double tp=(Tempreture(tmp.t,i)-Te)/mu[i];
	result.I=tmp.I+dt*(A[i]+B[i]*exp(-tp*tp))*tmp.I;
    result.lamb=1.0/result.I;
    return result;
}

void SolveOde(double T1_time, double End_time,int in){

    int i;
    if(solution[in]!=NULL) {
        free(solution[in]);
        solution[in]=NULL;
    }
    int n=ceil((End_time-T1_time)/dt)+1;
    solution[in]=malloc(sizeof(ODE)*n);
    solution[in][0].I=I0[in];
    solution[in][0].lamb=1.0/solution[0][in].I;
    solution[in][0].t=T1_time;
    for(i=1;i<n;i++){
        solution[in][i]=solver(solution[in][i-1],in);
    } 

  //  for(i=0;i<n;i++)
  //      printf("%e %e %e %e\n",solution[i].I,solution[i].t,solution[i].lamb);
}
