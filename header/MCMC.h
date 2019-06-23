#define MAXDATASIZE 2048
#define MAXMCTIMES 200000
#define INF 1e20
#define LOWBODUND -1e200
#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))

//global par

int seqlength=-1;
double lambda=0.001;  //mutation rate

#include "fasta.h"
#include "tree.h"
#include "randTree.h"
#include "init.h"
#include "ode.h"
#include "Coleselikelihood.h"
#include "Treelikelihood.h"

void MCMC(){

    int ncount=0;
    double ALPHA,r;
    double llhood;
    double pllhood=LOWBODUND;
    BiNode *TREE,*TEMPTREE;
    int id=getpid();
    char fname[30];

    double tl,tI0;
    double tmu,tA,tB;

    srand(time(0)*id);
    initTem();
    TREE=init();

    sprintf(fname,"result/%d_result.txt",id); // change here
    FILE *fp;
    fp=fopen(fname,"w+");

    printf("The initilization is over, then go to the calculation process!\n");
    
    while (ncount<MAXMCTIMES){

        tl=lambda;
        r=rand()/(RAND_MAX+0.0)-0.5;
        lambda*=exp(r);

        tmu=mu;
        r=rand()/(RAND_MAX+0.0)-0.5;
        mu*=exp(r*0.1);

        tA=A;
        r=rand()/(RAND_MAX+0.0)-0.5;
        A*=exp(r*0.5);

        tB=B;
        r=rand()/(RAND_MAX+0.0)-0.5;
        B*=exp(r*0.5);

        tI0=I0;
        r=rand()/(RAND_MAX+0.0)-0.5;
        I0*=exp(r*0.2);

        // random change tree
        TEMPTREE=cloneTree(TREE,NULL);
        if (ncount<0.2*MAXMCTIMES){
            RandTime(TEMPTREE);
            TEMPTREE=SPRtheTree(TEMPTREE);
        }
        else{
            RandTime(TEMPTREE);
            TEMPTREE=NNItheTree(TEMPTREE);
        }
        //calculate likelihood
        SolveOde(TEMPTREE->NodeTime,416);
        double llhood1=coallikelihood(TEMPTREE);
        double llhood2=LogTreeLikelihood(TEMPTREE);
        llhood=llhood1+llhood2;
        //printf("%lf %lf %d\n",tt,llhood-tt,ncount);
        ALPHA=MIN(0,llhood-pllhood);
        r=log(rand()/(RAND_MAX+0.0));

        if (r<ALPHA){
             pllhood=llhood;
             DestroyTree(TREE);
             TREE=cloneTree(TEMPTREE,NULL);
             DestroyTree(TEMPTREE);
             printf("%lf %lf %d\n",llhood1,llhood2,ncount);
            // if (ncount>10000) printf("%e\t%e\t%e\t%e\t%lf\t%lf\t%lf\t%lf\n",
            // lambda,A,B,C,I0,Te,mu,TREE->NodeTime);
             fprintf(fp,"%e\t%e\t%e\t%lf\t%lf\t%lf\n",
             lambda,A,B,I0,mu,TREE->NodeTime);
        }
        else{
            DestroyTree(TEMPTREE);
            lambda=tl;
            A=tA;
            B=tB;
            I0=tI0;
            mu=tmu;
        }
        
        ncount++;
       // printf("finished %dth llhood1=%lf 2=%lf\n",ncount,llhood1,llhood2);

    }
    fclose(fp);
    FILE *ft;
    sprintf(fname,"result/%d_Tree.txt",id);
    ft=fopen(fname,"w+");
    printTree(TREE,ft);
    fclose(ft);
    DestroyTree(TREE);

}

