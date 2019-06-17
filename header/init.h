typedef struct ComMat{
    double t;
    int nt;
}ComMat;
ComMat data[2];

double *TEM;

BiNode *init()
{

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int   L;
    char filename[100];
	strcpy(filename,"align_data/outcds.fasta");

	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	ffp = OpenFASTA(filename);

    if (ffp==NULL) {
        printf("wrong input!\n");
        return NULL;
    }
	data[0].t=416.0;
	data[0].nt=0;
	data[1].t=353.0;
	data[1].nt=0;

	while(ReadFASTA(ffp, &seq, &name, &L)){
        S[n]=calloc(1,sizeof(BiNode));
		S[n]->lchi=NULL;
		S[n]->rchi=NULL;
		S[n]->parent=NULL;
        S[n]->sequence=calloc(strlen(seq)+1,sizeof(char));
        S[n]->name=calloc(16,sizeof(char));
	    strcpy(S[n]->sequence,seq);
	    strcpy(S[n]->name,name);
		int tmp=atoi(name+9);
		if (tmp==2014){
			S[n]->NodeTime=416.0;
			data[0].nt++;
		}
		else if (tmp==2013){
			S[n]->NodeTime=353.0;
			data[1].nt++;
		}
		else{
			printf("No That Time %d !\n",tmp);
			exit(1);
		}
	
		//S[n]->NodeTime=0.0;
	    free(name);
	    free(seq);
		n++;
	}
	seqlength=strlen(S[1]->sequence);
	
	//printf("number of node=%d\n",n);
	top=n-1;
	//printf("length %d\n",seqlength);
	int i,j,len=n;
	int judge[MAXDATASIZE];
	//for (i=0;i<len;i++) printf("S[%d]=%d\n",i,S[i]->data);

	for (i=0;i<MAXDATASIZE;i++) judge[i]=i;
	//srand(time(0));
	int r;
	while(len>1){ //shuffle algorithm to build the Tree
		top++;
		S[top]=calloc(1,sizeof(BiNode));
		S[top]->name=calloc(10,sizeof(char));
		strcpy(S[top]->name,"AYNULL");
		S[top]->sequence=NULL;
		r=rand()%len;
		j=judge[r];
		S[top]->lchi=S[j];
        S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		//printf("j=%d\n",j);

		r=rand()%len;
		j=judge[r];
		S[top]->rchi=S[j];
		S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		judge[len]=top;
		len++;

		S[top]->NodeTime=MIN(S[top]->lchi->NodeTime,S[top]->rchi->NodeTime)-20;
		S[top]->parent=NULL;
		//printf("j=%d\n\n",j);

	}
    //printf("top=%d\n",top);
	//for (i=0;i<top+1;i++) printf("S[%d]=%d\n",i,S[i]->data);
	CloseFASTA(ffp);
	return S[top];
}

void initTem(){

	char filename[100];
	strcpy(filename,"align_data/weather_data.txt");
	FILE *fp;
	fp=fopen(filename,"r");
    char str[100];
    double a;
	int len=512,n=0;
	TEM=malloc(sizeof(double)*len);
    while(!feof(fp)){
        fscanf(fp,"%s   %lf\n",str,&a);
		TEM[n]=a;
		n++;
		if(n==len){
			len+=128;
			TEM=realloc(TEM,sizeof(double)*len);
		}
        
    }
	fclose(fp);
}
