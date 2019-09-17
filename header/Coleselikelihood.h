
int Comparator(const void * a, const void * b)
{
    return (*(ComMat *)a).t < (*(ComMat *)b).t ? 1 : -1;
}

double my_log(double x)
{
    double a=(1.0+x)/2.0;
    double b=sqrt(x);
    while(1){
        a=(a+b)/2.0;
        b=sqrt(a*b);
        if(fabs(a-b)<0.001) break;
    }
    return 2.0*(x-1.0)/(a+b);
}



double coallikelihood(BiNode *T,int in)
{
    
    ComMat *mlist;
    int len=256;
    mlist=(ComMat *)malloc(sizeof(ComMat)*len);
    int i,n=0;
    for(i=0;i<2;i++){
        mlist[n].t=data[i].t;
        mlist[n].nt=data[i].nt;
        n++;
    } 
    TimeListN=0;
    PrintTime(T);
    for(i=0;i<TimeListN;i++){
        if (n>=len){
            len+=256;
            mlist=(ComMat *)realloc(mlist,len*sizeof(ComMat));
        }
        mlist[n].t=TimeList[i];
        mlist[n].nt=-1;
        n++;
    }
    qsort(mlist,n,sizeof(ComMat),Comparator);
    double sums=0.0;
    double nnt=mlist[0].nt;
    double init_time=T->NodeTime;
    if (nnt<0) {
        printf("NUMBER WORONG!@ coleselikelihood.h lines 52\n");
        exit(1);
    }
    for(i=1;i<n;i++){
        if (nnt<2){
            if (mlist[i].nt==-1) nnt=2;
            else nnt=1;
        }
        int i1=floor((mlist[i].t-init_time)/dt);
        int i2=floor((mlist[i-1].t-init_time)/dt);
        int j;
        double lamv=0.0;
        for(j=i1;j<=i2;j++)
            lamv+=solution[in][j].lamb*dt;
        
        if (mlist[i].nt==-1){
            sums+=my_log(solution[in][i1].lamb*(nnt-1.0)*nnt/2.0)-lamv*(nnt-1.0)*nnt/2.0;
            // printf("1 %lf %lf %lf %lf\n",solution[i1].lamb1,sums,solution[i1].t,solution[i1].I1);
        }
        else{
            sums-=lamv*(nnt-1.0)*nnt/2.0;
        }
	    nnt+=mlist[i-1].nt;
    }
    free(mlist);
	return sums;
}
/*
double timelikelihood()
{
    return -extime*Nbegin*Ncodon*lambda*omega1;
}
*/
