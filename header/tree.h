typedef struct BiNode{
	struct BiNode *lchi,*rchi;
	struct BiNode *parent;
	char *name;
	char *sequence;
	double NodeTime;
}BiNode;

double TimeList[MAXDATASIZE]={0.0};
int TimeListN=0;

double logplus(double loga,double logb){

    if(loga>logb){
        return logb+log(1+exp(loga-logb));
    }else{
        return loga+log(1+exp(logb-loga));
    }

}

int printTree(BiNode* T,FILE *fp){
	if (T->lchi==NULL && T->rchi==NULL) {
		fprintf(fp,"%s",T->name);
		//printf(" %lf",T->NodeTime-T->parent->NodeTime);
		//printf("%c ",T->sequence[0]);
		return 1;
	}
	fprintf(fp,"(");
	//printf("%d(",T->data);
	printTree(T->lchi,fp);
   // printf(" %lf",T->lchi->NodeTime-T->NodeTime);
	fprintf(fp,",");
	printTree(T->rchi,fp);
   // printf(" %lf",T->lchi->NodeTime-T->NodeTime);
	fprintf(fp,")");
	return 1;
}
int PrintTime(BiNode *T){

    if (T->lchi==NULL && T->rchi==NULL){

         // TimeList[TimeListN]=T->NodeTime;
         // TimeListN++;
          return 1;

    }

    TimeList[TimeListN]=T->NodeTime;
    TimeListN++;
    PrintTime(T->lchi);
    PrintTime(T->rchi);
    return 1;
}

void DestroyTree(BiNode *T)
{
    if (T==NULL) return;

    DestroyTree(T->lchi);
    DestroyTree(T->rchi);

    if (T->sequence!=NULL) free(T->sequence);
    free(T->name);
    free(T);
}

int isTerminal(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    else
        return 0;
}
int isRoot(BiNode *T){

    if (T->parent==NULL)
        return 1;
    else
        return 0;
}


int TreeSize(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    if (T->lchi==NULL)
        return 1+TreeSize(T->rchi);
    if (T->rchi==NULL)
        return 1+TreeSize(T->lchi);
    return 1+TreeSize(T->lchi)+TreeSize(T->rchi);

}
int TreeLeaveSize(BiNode *T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    if (T->lchi==NULL)
        return TreeLeaveSize(T->rchi);
    if (T->rchi==NULL)
        return TreeLeaveSize(T->lchi);
    return TreeLeaveSize(T->lchi)+TreeLeaveSize(T->rchi);

}
BiNode *cloneTree(BiNode *T,BiNode *Tp){


    if (T==NULL) return NULL;

    BiNode *tmp;
    tmp=calloc(1,sizeof(BiNode));
    tmp->NodeTime=T->NodeTime;
    tmp->name=calloc(10,sizeof(char));
    strcpy(tmp->name,T->name);
    if(T->sequence==NULL){
        tmp->sequence=NULL;
    }
    else{
        tmp->sequence=calloc(strlen(T->sequence)+1,sizeof(char));
        strcpy(tmp->sequence,T->sequence);
    }
    tmp->parent=Tp;
    tmp->lchi=cloneTree(T->lchi,tmp);
    tmp->rchi=cloneTree(T->rchi,tmp);


    return tmp;

}

