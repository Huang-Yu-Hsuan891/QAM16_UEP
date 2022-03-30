#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

unsigned long long SEED =  312891; // have value //unsigned long long SEED = 3881111387891;
unsigned long long RANV;
int RANI = 0;
double Ranq1();
void normal(double sig, double *n1, double *n2);
double CHK(double L1, double L2);
double table[8] = {0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, 0};

double minimum(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8);
void table_symbol(int i, int a, int b, int c, int d);
void x_estimate_m(int i, double a, double b);
void table_receive(int i, double a, double b);

double **s;
int srow = 1844;  // Na / 4 transmitt how many symbols (transmitter)
int scolumn = 2; //1 symbol vector have 2 index

double **m_estimate;
int mrow = 1844;  // Na / 4 transmitt estimate symbols (receiver+estimate)
int mcolumn = 2; //1 symbol vector have 2 index

int **receive_sym;

int receivesymrow = 1844; // Na/4 transmitt how many symbols
int receivesymcolumn = 4; //16QAM PRESENT 1symbol = 4bits

int main() {
    // for declaration
    int i, j, k, m;                  //for counting
    int n, rc;                       // n is column = 1024 and rc is row = 512
    //int dv, dc;                      // dv: column have #1  = 3 ; and dc: row have #1 = 6
    int dv1, dv2, dc1, dc2;
    //int **L = NULL;          // check node connect 6 variable nodes
    ///int Llenrow = rc;        // rc = 512
    //int Llencolumn = dc;     // dc = 6
    
    int **L1 = NULL;                 // check node connect 6 variable nodes
    int L1lenrow = rc/2;        // rc = 3688/2
    int L1lencolumn = dc1;     // dc = 12
    int **L2 = NULL;                 // check node connect 6 variable nodes
    int L2lenrow = rc/2;        // rc = 3688/2
    int L2lencolumn = dc2;     // dc = 6

    //int **M = NULL;
    //int Mlenrow = n;       // n = 1024
    //int Mlencolumn = dv;   //dv = 3

    int **M1 = NULL;
    int M1lenrow = n/2;       // n = 7376/2
    int M1lencolumn = dv1;   // dv1 = 3
    int **M2 = NULL;
    int M2lenrow = n/2;       // n = 7376/2
    int M2lencolumn = dv2;   // dv2 = 6

    int error_block = 0;             // receive 100 error block
    int num = 0;           // how many blocks transmitted 
    double x, y;           // for noise
    
    int krc = 3688;                // krc = 512
    int **G;
    int Glenrow = krc;
    int Glencolumn = n;

    int *u;
    int ulen = krc;
    int temp;     // for changing the order od M and L
    double d1,d2,d3,d4,d5,d6,d7,d8;
    double s1,s2,s3,s4,s5,s6,s7,s8;


    

    // open file
    FILE *fpr;
    fpr=fopen("paritycheckmatrix.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d; row = %d\n", n, rc);
    fscanf(fpr,"%d",&dc1);
    fscanf(fpr,"%d",&dc2);
    fscanf(fpr,"%d",&dv1);
    fscanf(fpr,"%d",&dv2);
    
    

    L1lenrow = rc/2;
    L1lencolumn = dc1;
    L2lenrow = rc/2;
    L2lencolumn = dc2;

    M1lenrow = n/2; 
    M1lencolumn = dv1;
    M2lenrow = n/2; 
    M2lencolumn = dv2;

    M1 = (int **)malloc(M1lenrow * sizeof(int *));
    for (i = 0; i < M1lenrow; i++) M1[i] = (int *)malloc(M1lencolumn * sizeof(int));
    M2 = (int **)malloc(M1lenrow * sizeof(int *));
    for (i = 0; i < M2lenrow; i++) M2[i] = (int *)malloc(M2lencolumn * sizeof(int));
    L1 = (int **)malloc(L1lenrow * sizeof(int *));
    for (i = 0; i < L1lenrow; i++) L1[i] = (int *)malloc(L1lencolumn * sizeof(int));
    L2 = (int **)malloc(L2lenrow * sizeof(int *));
    for (i = 0; i < L2lenrow; i++) L2[i] = (int *)malloc(L2lencolumn * sizeof(int));
    printf("y\n");
printf("column = %d; row = %d\n", n, rc);
    for (j = 0; j < M1lenrow; j++) 
        for (i = 0; i < M1lencolumn; i++) 
            fscanf(fpr,"%d",&M1[j][i]);

    for (j = 0; j < M2lenrow; j++) 
        for (i = 0; i < M2lencolumn; i++) 
            fscanf(fpr,"%d",&M2[j][i]);

    for (j = 0; j < M1lenrow; j++)     // FOR ORDERING 
        for (i = 0; i < M1lencolumn; i++) 
            for (m = i; m < M1lencolumn; m++) 
                if (M1[j][m] < M1[j][i]) {
                    temp = M1[j][m];
                    M1[j][m] = M1[j][i];
                    M1[j][i] = temp;
                }
    for (j = 0; j < M2lenrow; j++) 
        for (i = 0; i < M2lencolumn; i++) 
            for (m = i; m < M2lencolumn; m++) 
                if (M2[j][m] < M2[j][i]) {
                    temp = M2[j][m];
                    M2[j][m] = M2[j][i];
                    M2[j][i] = temp;
                }

    for (i = 0; i < L1lenrow; i++)     // FOR ORDERING 
        for (j = 0; j < L1lencolumn; j++) 
            fscanf(fpr,"%d",&L1[i][j]);

    for (i = 0; i < L2lenrow; i++) 
        for (j = 0; j < L2lencolumn; j++) 
            fscanf(fpr,"%d",&L2[i][j]);
    for (i = 0; i < L1lenrow; i++) 
        for (j = 0; j < L1lencolumn; j++) 
            for (m = j; m < L1lencolumn; m++) 
                if (L1[i][m] < L1[i][j]) {
                    temp = L1[i][m];
                    L1[i][m] =L1[i][j];
                    L1[i][j] = temp;
                }

    for (i = 0; i < L2lenrow; i++) 
        for (j = 0; j < L2lencolumn; j++) 
            for (m = j; m < L2lencolumn; m++) 
                if (L2[i][m] < L2[i][j]) {
                    temp = L2[i][m];
                    L2[i][m] = L2[i][j];
                    L2[i][j] = temp;
                }

    
    fclose(fpr);
    // close file

    // OPEN FILE
    FILE *fpr1;
    fpr1=fopen("generatorC2.txt","r");
    printf("column = %d\n", n);
    printf("row = %d\n", krc);
    Glenrow = krc;
    Glencolumn = n;
    G = (int **)malloc(Glenrow * sizeof(int *));
    for (i = 0; i < Glenrow; i++) G[i] = (int *)malloc(Glencolumn * sizeof(int));
    for (i = 0; i < Glenrow; i++) 
        for (j = 0; j < Glencolumn; j++) 
            fscanf(fpr1,"%d",&G[i][j]);    
    fclose(fpr1);
     // close file

    ulen = krc; // ulen = 3688
    u = (int *)malloc(ulen * sizeof(int));
    
    
    int *codarray;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;
    codarraylen = n;
    codarray = (int *)malloc(codarraylen * sizeof(int));
    int **sym;         // present symbol ( 4bits to 1symbols)
    int symrow = 1844; // Na/4 transmitt how many symbols
    int symcolumn = 4; //16QAM PRESENT 1symbol = 4bits
    sym = (int **)malloc(symrow * sizeof(int *));
    for (i = 0; i < symrow; i++)
        sym[i] = (int *)malloc(symcolumn * sizeof(int));
    s = (double **)malloc(srow * sizeof(double *));
    for (i = 0; i < srow; i++)
        s[i] = (double *)malloc(scolumn * sizeof(double));
    
    double **x_receive;
    int xrow = 1844;  // Na / 4 receive how many symbols (receiver)
    int xcolumn = 2; //1 symbol vector have 2 index
    x_receive = (double **)malloc(xrow * sizeof(double *));
    for (i = 0; i < xrow; i++)
        x_receive[i] = (double *)malloc(xcolumn * sizeof(double));

    m_estimate = (double **)malloc(mrow * sizeof(double *));
    for (i = 0; i < mrow; i++)
        m_estimate[i] = (double *)malloc(mcolumn * sizeof(double));

    receive_sym = (int **)malloc(receivesymrow * sizeof(int *));
    for (i = 0; i < receivesymrow; i++)
        receive_sym[i] = (int *)malloc(receivesymcolumn * sizeof(int));
    int *recieve_codarray;
    recieve_codarray = (int *)malloc(codarraylen * sizeof(int));


    double ebn0s = 0.2;//SNR
    double ebn0;
    double sigma;
    //sigma = sqrt(5.0 / (pow(10, ebn0s / 10)));
    sigma = sqrt(1.25 * 2 /* 1/R = 1/0.5 */ / (pow(10, ebn0s / 10)));
    ebn0 = pow(10, ebn0s/10);
    printf("%g,%g \n", pow(10, ebn0s / 10), sigma);

    // for CODE part
    // receive 100 error block
    double *Lj;             // LLR; LENGTH = 7376
    int Ljlen = n;          // LENGTH = 7376
    Lj = (double *)malloc(Ljlen * sizeof(double));
    double **qij1 = NULL;    // from down to top
    int qij1row = dv1;        // qijrow = dv1 =3
    int qij1column = n/2;      // qijcolumn = n = 7376/2
    double **qij2 = NULL;    // from down to top
    int qij2row = dv2;        // qijrow = dv = 6
    int qij2column = n/2;      // qijcolumn = n = 7376/2
    double **uij1 = NULL;    // from top to down
    int uij1row = rc/2;        // uijrow = rc = 3688/2
    int uij1column = dc1;     // uijcloumn = dc 12
    double **uij2 = NULL;    // from top to down
    int uij2row = rc/2;        // uijrow = rc = 3688/2
    int uij2column = dc2;     // uijcloumn = dc 6
    
    qij1row = dv1;        // qijrow = dv = 3
    qij1column = n/2;      // qijcolumn = n = 7376/2
    qij1 = (double **)malloc(qij1row * sizeof(double *));
    for (i = 0; i < qij1row; i++) qij1[i] = (double *)malloc(qij1column * sizeof(double));
    qij2row = dv2;        // qijrow = dv = 3
    qij2column = n/2;      // qijcolumn = n = 7376/2
    qij2 = (double **)malloc(qij2row * sizeof(double *));
    for (i = 0; i < qij2row; i++) qij2[i] = (double *)malloc(qij2column * sizeof(double));

    uij1row = rc/2;        // uijrow = rc = 512
    uij1column = dc1;     // uijcolumn = dc = 6 
    uij1 = (double **)malloc(uij1row * sizeof(double *));
    for (i = 0; i < uij1row; i++) uij1[i] = (double *)malloc(uij1column * sizeof(double));
    uij2row = rc/2;        // uijrow = rc = 512
    uij2column = dc1;     // uijcolumn = dc = 6 
    uij2 = (double **)malloc(uij2row * sizeof(double *));
    for (i = 0; i < uij2row; i++) uij2[i] = (double *)malloc(uij2column * sizeof(double));
    
    int restart = 0;
    double tempqij1[11]; 
    double tempqij2[5]; 
    double tempuij;
    double temp1uij1[3];
    double temp1uij2[6];
    double temp1qij; 
    int valL;
    int valL2;
    int *comput;
    int computlen = n;      // computlen = 1024
    int *comput1;
    int comput1len = rc;    // comput1len = 512
    int *comput2;
    int comput2len = rc;
    computlen = n;
    comput = (int *)malloc(computlen * sizeof(int));
    comput1len = rc;
    comput1 = (int *)malloc(comput1len * sizeof(int));
    comput2len = rc;
    comput2 = (int *)malloc(comput2len * sizeof(int));
    double *qj;
    int qjlen = n;     //qjlen = 1024
    int *checkbit;
    int checkbitlen = rc;   // checkbitlen = 512
    qjlen = n;
    qj = (double *)malloc(qjlen * sizeof(double));
    checkbitlen = rc;
    checkbit = (int *)malloc(checkbitlen * sizeof(int));
    int *output;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;
    outputarray = n;
    output = (int *)malloc( outputarray * sizeof(int) );
    int error1 = 0;
    int totalerror1 = 0;
    int error2 = 0;
    int totalerror2 = 0;
    double BER1;

    double BER2;
    double FER1;
    double FER2;

    while (error_block < 100/*100*/) {
        for (i = 0; i < codarraylen; i++) codarray[i] = 0;
        //printf("yes");
        if (num == 0) {
            u[0] = 1;
            u[1] = 0;
            u[2] = 0;
            u[3] = 0;
            u[4] = 0;
            u[5] = 0;
            for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
        } else {
            u[0] = (u[krc-5] + u[krc-6]) % 2;
            u[1] = (u[krc-4] + u[krc-5]) % 2;
            u[2] = (u[krc-3] + u[krc-4]) % 2;
            u[3] = (u[krc-2] + u[krc-3]) % 2;
            u[4] = (u[krc-1] + u[krc-2]) % 2;
            u[5] = (u[krc-0] + u[krc-1]) % 2;
            for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
        }//for (i = 0; i < krc; i++) printf("u[%d] = %d ;", i, u[i]);
        num++;
        printf("codarraylen = %d ",codarraylen);
        for (i = 0; i < ulen; i++) 
            if (u[i] == 1) {
                for (j = 0; j < n; j++) codarray[j] = (codarray[j] + G[i][j]) % 2;
            }//for (i = 0; i < codarraylen; i++) printf("codarray[%d] = %d ;", i, codarray[i]);
        for (i = 0; i < symrow; i++)
            for (j = 0; j < symcolumn; j++) {
                sym[i][j] = codarray[4 * i + j];
                //if (i == 0 || i == 10) printf("sym[%d][%d] = %d; ", i, j, sym[i][j]);
            }
        for (i = 0; i < srow; i++)
            table_symbol(i, sym[i][0], sym[i][1], sym[i][2], sym[i][3]);
        for(i = 0; i < srow; i++) {
            normal(sigma, &x, &y);
            x_receive[i][0] = s[i][0] + x;
            x_receive[i][1] = s[i][1] + y;
        }
        /*for (i = 0; i < xrow; i++) x_estimate_m(i,x_receive[i][0],x_receive[i][1]);
        for (i = 0; i < receivesymrow; i++) {
            table_receive(i, m_estimate[i][0],m_estimate[i][1]);
            if ( i == 0 || i == 1) printf("receive_sym[%d]= %d %d %d %d\n",i,receive_sym[i][0],receive_sym[i][1],receive_sym[i][2],receive_sym[i][3]);
        }
        for (i = 0; i < codarraylen; i++) recieve_codarray[i] = receive_sym[i/4][i%4];*/
        
        for (i = 0; i < Ljlen; i++) {
            //printf("i = %d\n",i);
            if ((i % 4) == 0) {
                //FOR L1
                d1 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 1),2);         // X1 = 0
                d2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 1),2);
                d3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 1),2);
                d4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 1),2);
                d5 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 3),2);
                d6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 3),2);
                d7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 3),2);
                d8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 3),2);
                s1 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 1),2);         // X1 = 1
                s2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 1),2);
                s3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 1),2);
                s4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 1),2);
                s5 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 3),2);
                s6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 3),2);
                s7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 3),2);
                s8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 3),2);
                //printf("%g,%g,%g,%g,%g,%g,%g,%g\n",d1,d2,d3,d4,d5,d6,d7,d8);
                Lj[i] = 0.4 * 0.5 * ebn0 * (- minimum(d1,d2,d3,d4,d5,d6,d7,d8) + minimum(s1,s2,s3,s4,s5,s6,s7,s8));    // L1 : for first bit in symbol 
                //printf("%g %g \n",minimum(d1,d2,d3,d4,d5,d6,d7,d8),minimum(s1,s2,s3,s4,s5,s6,s7,s8));
            } else if ((i % 4) == 1) {
                //FOR L2
                d1 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 1),2);         // X2 = 0
                d2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 1),2);
                d3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 1),2);
                d4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 1),2);
                d5 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 1),2);
                d6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 1),2);
                d7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 1),2);
                d8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 1),2);

                s1 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 3),2);         // X2 = 1
                s2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 3),2);
                s3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 3),2);
                s4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 3),2);
                s5 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 3),2);
                s6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 3),2);
                s7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 3),2);
                s8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 3),2);
                Lj[i] = 0.4 * 0.5 * ebn0 * (- minimum(d1,d2,d3,d4,d5,d6,d7,d8) + minimum(s1,s2,s3,s4,s5,s6,s7,s8));    // L2 : for first bit in symbol 
            } else if ((i % 4) == 2) {
                //FOR L3
                d1 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 1),2);         // X3 = 0
                d2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 1),2);
                d3 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 1),2);
                d4 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 1),2);
                d5 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 3),2);
                d6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 3),2);
                d7 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 3),2);
                d8 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 3),2);

                s1 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 1),2);         // X3 = 1
                s2 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 1),2);
                s3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 1),2);
                s4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 1),2);
                s5 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 3),2);
                s6 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 3),2);
                s7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 3),2);
                s8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 3),2);
                Lj[i] = 0.4 * 0.5 * ebn0 * (- minimum(d1,d2,d3,d4,d5,d6,d7,d8) + minimum(s1,s2,s3,s4,s5,s6,s7,s8));    // L3 : for first bit in symbol 
            } else if ((i % 4) == 3) {
                //FOR L4
                d1 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 1),2);         // X4 = 0
                d2 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 1),2);
                d3 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 1),2);
                d4 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 1),2);
                d5 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] + 3),2);
                d6 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] + 3),2);
                d7 = pow((x_receive[i / 4][0] - 1),2) + pow((x_receive[i/4][1] - 3),2);
                d8 = pow((x_receive[i / 4][0] + 1),2) + pow((x_receive[i/4][1] - 3),2);

                s1 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 1),2);         // X4 = 1
                s2 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 1),2);
                s3 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 1),2);
                s4 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 1),2);
                s5 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] + 3),2);
                s6 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] + 3),2);
                s7 = pow((x_receive[i / 4][0] + 3),2) + pow((x_receive[i/4][1] - 3),2);
                s8 = pow((x_receive[i / 4][0] - 3),2) + pow((x_receive[i/4][1] - 3),2);
                Lj[i] = 0.4 * 0.5 * ebn0 * (- minimum(d1,d2,d3,d4,d5,d6,d7,d8) + minimum(s1,s2,s3,s4,s5,s6,s7,s8));    // L4 : for first bit in symbol 
            } else Lj[i] = 0;
        } // END FOR LLR

        // the interative decoding algotrithm
        for (j = 0; j < qij1column; j++)                             // initialization
            for (i = 0; i < qij1row; i++) qij1[i][j] = Lj[j];
        for (j = 0; j < qij2column; j++)                             // initialization
            for (i = 0; i < qij2row; i++) qij2[i][j] = Lj[3688 + j];

        // message passing, for predetermined threshold = 200
        for (k = 0; k < 200 && restart != rc; k++) {
            restart = 0;
            //printf("k = %d ",k);
            // bottom-up
            for (i = 0; i < 11; i++) tempqij1[i] = 0.0;
            for (i = 0; i < 5; i++) tempqij2[i] = 0.0;
            for (i = 0; i < computlen; i++) comput[i] = 0;
            for (i = 0; i < rc; i++) {
                    if (i < (rc/2)) {
                        for (j = 0; j < 12; j++) {
                            for (m = 0; m < 11; m++) {
                                if (m < j) {
                                    valL = L1[i][m]-1;
                                    if (valL < (n/2)) tempqij1[m] = qij1[comput[valL]][valL];
                                    else tempqij1[m] = qij2[comput[valL]][valL - 3688];
                                    //if (i == 0&& j==8)printf("valL = %d ; ", valL);
                                } 
                                else if (m >= j) {
                                    valL = L1[i][m+1]-1;
                                    if (valL < (n/2)) tempqij1[m] = qij1[comput[valL]][valL];
                                    else tempqij1[m] = qij2[comput[valL]][valL - 3688];
                                    //if (i == 0&& j==8)printf("valL = %d ; ", valL);
                                }
                            }
                            tempuij = tempqij1[0];
                            //if (i == 0&& j==8) printf("tempuij = %.8g; ",tempuij);
                            for(m = 1; m < 11; m++) {
                                tempuij = CHK(tempuij, tempqij1[m]);
                                //if (i == 0&& j==8) printf("tempuij = %.8g; ",tempuij);
                                //if (i == 0&& j==8&&m==10) printf("m=10,tempqij1[m]=%.8g\n ",tempqij1[m]);
                            }
                            uij1[i][j] = tempuij;
                            //if (i == 0 && j==8) {
                                //printf("uij1[%d][%d]= %.7g ;",i,j,uij1[i][j]);
                                //for (m=0; m <11;m++)printf("tempqij1[%d] = %g; ",m,tempqij1[m]);
                            //}
                            //if (i == 1 && j == 8) printf("uij1[%d][%d] = %g;    ",i,j,uij1[i][j]);
                        }
                    } else {
                        for (j = 0; j < 6; j++) {
                            for (m = 0; m < 5; m++) {
                                if (m < j) {
                                    valL = L2[i-1844][m]-1;
                                    if (valL < (n/2)) tempqij2[m] = qij1[comput[valL]][valL];
                                    else tempqij2[m] = qij2[comput[valL]][valL-3688];
                                } 
                                else if (m >= j) {
                                    valL = L2[i-1844][m+1]-1;
                                    if (valL < (n/2)) tempqij2[m] = qij1[comput[valL]][valL];
                                    else tempqij2[m] = qij2[comput[valL]][valL-3688];
                                }
                            }
                            tempuij = tempqij2[0];
                            for(m = 1; m < 5; m++) {
                                tempuij = CHK(tempuij, tempqij2[m]);
                            }
                            uij2[i-1844][j] = tempuij;
                        }
                    }
                    if (i < (rc/2)) {
                        for (m = 0; m < 12; m++) {
                            comput[L1[i][m] - 1] += 1;
                        }
                    } else {
                        for (m = 0; m < 6; m++) {
                            comput[L2[i-1844][m] - 1] += 1;
                        }
                    }
                }

                // top-down
                for(i = 0; i < 3; i++) {
                    temp1uij1[i] = 0.0;
                }
                for(i = 0; i < 6; i++) {
                    temp1uij2[i] = 0.0;
                }
                for (i = 0; i < comput1len; i++) comput1[i] = 0;
                for (j = 0; j < n; j++) {
                    if (j < (n/2)) {
                        for (i = 0; i < 3; i++) {
                            for (m = 0; m < 2; m++) {
                                if (m < i) { 
                                    valL = M1[j][m] - 1;
                                    if (valL < (rc/2)) temp1uij1[m] = uij1[valL][comput1[valL]];
                                    else  temp1uij1[m] = uij2[valL-1844][comput1[valL]];
                                }
                                else if (m >= i) {
                                    valL = M1[j][m + 1] - 1;
                                    if (valL < (rc/2)) temp1uij1[m] = uij1[valL][comput1[valL]];
                                    else  temp1uij1[m] = uij2[valL-1844][comput1[valL]];
                                }
                            }
                            temp1uij1[2] = Lj[j];
                            qij1[i][j] = temp1uij1[0] + temp1uij1[1] + temp1uij1[2];
                            //if(/*comput1[valL] == 0&&*/(j == 901 || j == 1363 ||j == 1803 || j== 2664||j==3064||j==3423||j==4378||j==4872||j==5101||j==6315||j==6607||j==6929))printf("temp1uij1[0] = %g; temp1uij1[1] = %g; temp1uij1[2] = %g; \n",temp1uij1[0],temp1uij1[1],temp1uij1[2]);
                            //if(/*comput1[valL] == 0&&*/(j == 901 || j == 1363 ||j == 1803 || j== 2664||j==3064||j==3423||j==4378||j==4872||j==5101||j==6315||j==6607||j==6929)) printf("qij1[%d][%d] = %g ;\n",i,j,qij1[i][j]);
                        }
                    } else {
                        for (i = 0; i < 6; i++) {
                            for (m = 0; m < 5; m++) {
                                if (m < i) { 
                                    valL = M2[j-3688][m] - 1;
                                    if (valL < (rc/2)) temp1uij2[m] = uij1[valL][comput1[valL]]; 
                                    else temp1uij2[m] = uij2[valL-1844][comput1[valL]]; 
                                }
                                else if (m >= i) {
                                    valL = M2[j-3688][m + 1] - 1;
                                    if (valL < (rc/2)) temp1uij2[m] = uij1[valL][comput1[valL]]; 
                                    else temp1uij2[m] = uij2[valL-1844][comput1[valL]];
                                }
                            }
                            temp1uij2[5] = Lj[j];
                            qij2[i][j-3688] = temp1uij2[0] + temp1uij2[1] + temp1uij2[2] +temp1uij2[3]+temp1uij2[4]+temp1uij2[5];
                            //if((j == 901 || j == 1363 ||j == 1803 || j== 2664||j==3064||j==3423||j==4378||j==4872||j==5101||j==6315||j==6607||j==6929))printf("temp1uij2[0] = %g; temp1uij2[1] = %g; temp1uij2[2] = %g; temp1uij2[3] = %g;temp1uij2[4] = %g;temp1uij2[5] = %g;\n",temp1uij2[0],temp1uij2[1],temp1uij2[2],temp1uij2[3],temp1uij2[4],temp1uij2[5]);
                            //if((j == 901 || j == 1363 ||j == 1803 || j== 2664||j==3064||j==3423||j==4378||j==4872||j==5101||j==6315||j==6607||j==6929)) printf("qij2[%d][%d] = %g ;\n",i,j,qij2[i][j-3688]);
                        }
                    }
                    if (j < n/2) {
                        for (m = 0; m < 3; m++) {
                            comput1[M1[j][m] - 1] += 1;
                        }
                    } else {
                        for (m =0; m < 6; m++) {
                            comput1[M2[j-3688][m] - 1] += 1;
                        }
                    }
                }

                // decision
                for (i = 0; i < comput2len; i++) comput2[i] = 0;
                for (j = 0; j < n; j++) {
                    qj[j] = Lj[j];
                    if (j < (n/2)) {
                        for (i = 0; i < 3; i++) {
                            valL = M1[j][i] - 1;
                            if (valL < (rc/2)) qj[j] += uij1[valL][comput2[valL]];
                            else qj[j] += uij2[valL-1844][comput2[valL]];
                        }
                        if (qj[j] >= 0) output[j] = 0;
                        else if (qj[j] < 0) output[j] = 1;
                        for (i = 0; i < 3; i++) {
                            comput2[M1[j][i] - 1] += 1;
                        }
                    } else {
                        for (i = 0; i < 6; i++) {
                            valL = M2[j-3688][i] - 1;
                            if (valL < (rc/2)) qj[j] += uij1[valL][comput2[valL]];
                            else qj[j] += uij2[valL-1844][comput2[valL]];
                            //qj[j] += uij2[valL][comput2[valL]];
                        }
                        if (qj[j] >= 0) output[j] = 0;
                        else if (qj[j] < 0) output[j] = 1;
                        for (i = 0; i < 6; i++) {
                            comput2[M2[j-3688][i] - 1] += 1;
                        }
                    }
                    //if (j ==0 || j == 515|| j==4789||j ==7000) printf("qj[%d] = %g\n",j,qj[j]);
                }

                // to check Hx=0     
                for (i = 0; i < rc; i++) {
                    checkbit[i] = 0;
                    if (i < (rc/2)) {
                        for (j = 0; j < 12; j++) {
                            checkbit[i] += output[L1[i][j] - 1];
                        }
                        checkbit[i] = checkbit[i] % 2;
                    } else {
                        for (j = 0; j < 6; j++) {
                            checkbit[i] += output[L2[i-1844][j] - 1];
                        }
                        checkbit[i] = checkbit[i] % 2;
                    }
                }

                for (i = 0; i < rc; i++) {
                    if (checkbit[i] == 0) restart += 1; // restart = 408 is success
                }
                int restart1 = 0;
                int restart2 = 0;
                for (i = 0; i < rc/2; i++) {
                    if (checkbit[i] == 0) restart1 += 1; // restart = 408 is success
                }
                for (i = rc/2; i < rc; i++) {
                    if (checkbit[i] == 0) restart2 += 1; // restart = 408 is success
                }
            //printf("restart =%d\n",restart);
        
        }// end for iteration
        printf("errorblock = %d num %d \n",error_block,num);
        error1 = 0;
        error2 = 0;
        for(i = 0; i < n/2; i++) if (output[i] != codarray[i]) error1 += 1;
        for(i = n/2; i < n; i++) if (output[i] != codarray[i]) error2 += 1;
        if (error1 != 0) {
            error_block++;
            printf("error1 = %d \n",error1);
        }
        restart = 0;
        totalerror1 += error1;
        totalerror2 += error2;
    }
    BER1 = (double)totalerror1/(num*n);
    BER2 = (double)totalerror2/(num*n);
    FER1 = (double)100.0/num;
    //FER2 = (double)100.0/num;
    printf("totalerror1 = %d, BER1 = %g, FER = %g\n",totalerror1,BER1,FER1);
    printf("totalerror2 = %d, BER2 = %g, FER = %g\n",totalerror2,BER2,FER1);
    FILE *outfp1; 
        outfp1 = fopen("QAM16_UEP_LEVEL1.txt","a");
        fprintf(outfp1,"SNR = %g ; BER1 = %g ;BER2 = %g ; FER = %g\n", ebn0s, BER1,BER2,FER1);
    fclose(outfp1);
    return 0;
}

double Ranq1(){
    if (RANI == 0)
    {
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

void normal(double sig, double *n1, double *n2)
{
    double x1, x2;
    double s;
    //printf("sigma = %g\n", sig);
    do
    {
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s)) / s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s)) / s);
}

void table_symbol(int i, int a, int b, int c, int d)
{
    if (a == 1 && b == 1 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 0 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = 1;
    }
    else if (a == 0 && b == 0 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 1 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = -3;
    }
    else
    {
        s[i][0] = 0;
        s[i][1] = -0;
    }
}
void x_estimate_m(int i, double a, double b){
    if (a >= 2) m_estimate[i][0] = 3;
    else if (a < 2 && a >= 0) m_estimate[i][0] = 1;
    else if (a < 0 && a >= -2) m_estimate[i][0] = -1;
    else if (a < -2) m_estimate[i][0] = -3;
    else m_estimate[i][0] = 0;
    if (b >= 2) m_estimate[i][1] = 3;
    else if (b < 2 && b >= 0) m_estimate[i][1] = 1;
    else if (b < 0 && b >= -2) m_estimate[i][1] = -1;
    else if (b < -2) m_estimate[i][1] = -3;
    else m_estimate[i][1] = 0;
}
void table_receive(int i, double a, double b) {
    if (a == -3 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0; 
    }
}

double minimum(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
    double min = a1;
    if (min > a2) min = a2;
    if (min > a3) min = a3;
    if (min > a4) min = a4;
    if (min > a5) min = a5;
    if (min > a6) min = a6;
    if (min > a7) min = a7;
    if (min > a8) min = a8;
    return min;
}
double CHK(double L1, double L2) 
{   
    double sgn1, sgn2, min;

    if(L1>0) sgn1 = 1;
    else if(L1 == 0) sgn1 = 0;
    else sgn1 = -1;
    if(L2>0) sgn2 = 1;
    else if(L2 == 0) sgn2 = 0;
    else sgn2 = -1;
    if(fabs(L1) >= fabs(L2)) min = fabs(L2);
    else min = fabs(L1);
    double temp1,temp2;
    double ope1,ope2;
    temp1 = L1+L2;
    temp2 = L1-L2;
    ope1 = fabs(temp1);
    ope2 = fabs(temp2);
    double lnup, lndown;
    if (ope1 >= 0 && ope1 < 0.196) lnup = table[0];
    else if (ope1 >= 0.196 && ope1 < 0.433) lnup = table[1];
    else if (ope1 >= 0.433 && ope1 < 0.71) lnup = table[2];
    else if (ope1 >= 0.71 && ope1 < 1.05) lnup = table[3];
    else if (ope1 >= 1.05 && ope1 < 1.508) lnup = table[4];
    else if (ope1 >= 1.508 && ope1 < 2.252) lnup = table[5];
    else if (ope1 >= 2.252 && ope1 < 4.5) lnup = table[6];
    else if (ope1 >= 4.5) lnup = table[7];

    if (ope2 >= 0 && ope2 < 0.196) lndown = table[0];
    else if (ope2 >= 0.196 && ope2 < 0.433) lndown = table[1];
    else if (ope2 >= 0.433 && ope2 < 0.71) lndown = table[2];
    else if (ope2 >= 0.71 && ope2 < 1.05) lndown = table[3];
    else if (ope2 >= 1.05 && ope2 < 1.508) lndown = table[4];
    else if (ope2 >= 1.508 && ope2 < 2.252) lndown = table[5];
    else if (ope2 >= 2.252 && ope2 < 4.5) lndown = table[6];
    else if (ope2 >= 4.5) lndown = table[7];

    return sgn1 * sgn2 * min +lnup-lndown /*answer*/;
}