#define MAXDATASIZE 2048
#define MAXMCTIMES 600000
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
    BiNode *TREE[2],*TEMPTREE[2];
    int id=getpid();
    char fname[30];

    double tl,tI0[2],tTe;
    double tmu[2],tA[2],tB[2];

    srand(time(0)*id);
    initTem0();
    initTem1();

    while(1){
        TREE[0]=init(0);
        TREE[1]=init(1);
        if (TREE[0]->NodeTime>100 && TREE[1]->NodeTime>100) break;
        else{
            DestroyTree(TREE[0]);
            DestroyTree(TREE[1]);
        }
    }
    sprintf(fname,"result/%d_result.txt",id); // change here
    FILE *fp;
    fp=fopen(fname,"w+");

    printf("The initilization is over, then go to the calculation process!\n");
    
    while (ncount<MAXMCTIMES){

        tl=lambda;
        r=rand()/(RAND_MAX+0.0)-0.5;
        lambda*=exp(r);

        tTe=Te;
        r=rand()/(RAND_MAX+0.0)-0.5;
        Te*=exp(r*0.1);

        int xt;
        for(xt=0;xt<2;xt++){
            tA[xt]=A[xt];
            r=rand()/(RAND_MAX+0.0)-0.5;
            A[xt]*=exp(r*0.5);

            tB[xt]=B[xt];
            r=rand()/(RAND_MAX+0.0)-0.5;
            B[xt]*=exp(r*0.5);

            tI0[xt]=I0[xt];
            r=rand()/(RAND_MAX+0.0)-0.5;
            I0[xt]*=exp(r*0.3);

	    tmu[xt]=mu[xt];
            r=rand()/(RAND_MAX+0.0)-0.5;
            mu[xt]*=exp(r*0.1);
        }

        // random change tree
       
        for(xt=0;xt<2;xt++){
            TEMPTREE[xt]=cloneTree(TREE[xt],NULL);
            if (ncount<0.2*MAXMCTIMES){
                RandTime(TEMPTREE[xt]);
                TEMPTREE[xt]=SPRtheTree(TEMPTREE[xt]);
            }
            else{
                RandTime(TEMPTREE[xt]);
                TEMPTREE[xt]=NNItheTree(TEMPTREE[xt]);
            }
            //calculate likelihood
            SolveOde(TEMPTREE[xt]->NodeTime,416,xt);
        }

        double llhood1=coallikelihood(TEMPTREE[0],0)+coallikelihood(TEMPTREE[1],1);
        double llhood2=LogTreeLikelihood(TEMPTREE[0])+LogTreeLikelihood(TEMPTREE[1]);
        llhood=llhood1+llhood2;
        //printf("%lf %lf %d\n",tt,llhood-tt,ncount);
        ALPHA=MIN(0,llhood-pllhood);
        r=log(rand()/(RAND_MAX+0.0));

        if (r<ALPHA){
            pllhood=llhood;
            for(xt=0;xt<2;xt++){
                DestroyTree(TREE[xt]);
                TREE[xt]=cloneTree(TEMPTREE[xt],NULL);
                DestroyTree(TEMPTREE[xt]);
            }
            printf("%lf %lf %d\n",llhood1,llhood2,ncount);
            // if (ncount>10000) printf("%e\t%e\t%e\t%e\t%lf\t%lf\t%lf\t%lf\n",
            // lambda,A,B,C,I0,Te,mu,TREE->NodeTime);
            fprintf(fp,"%e\t%lf",lambda,Te);
            for(xt=0;xt<2;xt++){
                 fprintf(fp,"\t%e\t%e\t%e\t%lf\t%lf",A[xt],B[xt],I0[xt],mu[xt],TREE[xt]->NodeTime);
            }
            fprintf(fp,"\n");
        }
        else{
            for(xt=0;xt<2;xt++){
                DestroyTree(TEMPTREE[xt]);
                A[xt]=tA[xt];
                B[xt]=tB[xt];
                I0[xt]=tI0[xt];         
		mu[xt]=tmu[xt];
            }
            lambda=tl;
            Te=tTe;
        }
        
        ncount++;
       // printf("finished %dth llhood1=%lf 2=%lf\n",ncount,llhood1,llhood2);

    }
    fclose(fp);
    FILE *ft;
    sprintf(fname,"result/%d_Tree.txt",id);
    ft=fopen(fname,"w+");
    int j;
    for(j=0;j<2;j++){
        printTree(TREE[j],ft);
        fprintf(ft,"\n\n");
        DestroyTree(TREE[j]);
    }
    fclose(ft);

}

