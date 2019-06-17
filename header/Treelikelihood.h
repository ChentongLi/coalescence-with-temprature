typedef struct LNode{
    
    double TCAG[4];
    
}LNode;

double ProbJC69(int i, int j, double t){
    if (i==j) return 0.25+0.75*exp(-4.0*lambda*t);
    else  return 0.25-0.25*exp(-4.0*lambda*t);
}


LNode TreeLikelihood(BiNode *T, int n){ //the time here is the molecular clock
    
    LNode like;
    
    like.TCAG[0]=0.0; like.TCAG[1]=0.0; like.TCAG[2]=0.0; like.TCAG[3]=0.0;
    
    if (isTerminal(T)){
        if (T->sequence[n]=='T') {like.TCAG[0]=1.0; return like;}
        if (T->sequence[n]=='C') {like.TCAG[1]=1.0; return like;}
        if (T->sequence[n]=='A') {like.TCAG[2]=1.0; return like;}
        if (T->sequence[n]=='G') {like.TCAG[3]=1.0; return like;}
        if (T->sequence[n]=='R') {like.TCAG[2]=0.5; like.TCAG[3]=0.5; return like;}
        if (T->sequence[n]=='Y') {like.TCAG[0]=0.5; like.TCAG[1]=0.5; return like;}
        if (T->sequence[n]=='M') {like.TCAG[1]=0.5; like.TCAG[2]=0.5; return like;}
        if (T->sequence[n]=='K') {like.TCAG[0]=0.5; like.TCAG[3]=0.5; return like;}
        printf("NO NUCLITIDE WRONG! at the place %d and this is %c at tree 1 named %s\n",n,T->sequence[n],T->name);
        //printf("%s\n\n",T->sequence);
        printf("%s\n",T->sequence+n);
        exit(1);
    }
    
    LNode plikeLeft,plikeRight;
    int i,j;
    double tLeft=T->lchi->NodeTime - T->NodeTime;
    double tRight=T->rchi->NodeTime - T->NodeTime;
    
    if (tLeft<0 || tRight<0) {
        printf("time wrong!\n");
        exit(1);
    }
    
    plikeLeft=TreeLikelihood(T->lchi,n);
    plikeRight=TreeLikelihood(T->rchi,n);
    
    for (i=0;i<4;i++)
        for (j=0;j<4;j++)
            like.TCAG[i]+=ProbJC69(i,j,tLeft)*plikeLeft.TCAG[j];
    double temp=0.0;
    
    for (i=0;i<4;i++){
        for (j=0;j<4;j++)
            temp+=ProbJC69(i,j,tRight)*plikeRight.TCAG[j];
        
        like.TCAG[i]=like.TCAG[i]*temp;
        //temp=0.0;
    }
    
    return like;
}


double LogTreeLikelihood(BiNode *root){// p(D|G,t)
    
    double m=0.0,m2;
    int i,j;
    int n=0;
    LNode temp;
    //printf("%d\n",seqlength);
    for(i=0;i<seqlength;i++){
        temp=TreeLikelihood(root,i);
        m2=temp.TCAG[0]+temp.TCAG[1]+temp.TCAG[2]+temp.TCAG[3];
        m2=m2*0.25;
        //printf("m2=%lf\n",m2);
        m+=log(m2);
    }
    return m;
    
}
