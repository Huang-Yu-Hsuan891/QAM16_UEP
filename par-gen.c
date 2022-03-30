#include <stdio.h>
#include <stdlib.h>

int main() {
    // for H1
    int n, rc;
    int e;
    int n1, rc1;
    // for H2
    int N, RC;
    int E;
    int N1, RC1;
    int i, j, m, k;
    int **H1Cmask;
    int H1Cmaskrow, H1Cmaskcolumn;
    int **H2Cmask;
    int H2Cmaskrow, H2Cmaskcolumn;
    int **H1C;
    int H1Crow, H1Ccolumn;
    int **H2C;
    int H2Crow, H2Ccolumn;


    int **H;
    int Hrow, Hcolumn;
    int **H1;
    int H1row, H1column;


    FILE *fpr;
    fpr = fopen("H1.txt", "r");
    fscanf(fpr,"%d",&e);
    fscanf(fpr, "%d",&n);
    fscanf(fpr, "%d",&rc);
    fscanf(fpr, "%d",&rc1);
    fscanf(fpr, "%d",&n1);
    n1 = n1 * 2;
    H1Cmaskrow = rc1;
    H1Cmaskcolumn = n1;
    H1Cmask = (int **)malloc(H1Cmaskrow * sizeof(int *));
    for (i = 0; i < H1Cmaskrow; i++) H1Cmask[i] = (int *)malloc(H1Cmaskcolumn * sizeof(int));
    for (i = 0; i < rc1; i++) {
        for (j = 0; j < n1; j++) {
            fscanf(fpr, "%d", &H1Cmask[i][j]);
        }
    }
    fclose(fpr);

    FILE *fpr2;
    fpr2 = fopen("H2.txt", "r");
    fscanf(fpr,"%d",&E);
    fscanf(fpr, "%d",&N);
    fscanf(fpr, "%d",&RC);
    fscanf(fpr, "%d",&RC1);
    fscanf(fpr, "%d",&N1);
    N1 = N1 * 2;
    H2Cmaskrow = RC1;
    H2Cmaskcolumn = N1;
    H2Cmask = (int **)malloc(H2Cmaskrow * sizeof(int *));
    for (i = 0; i < H2Cmaskrow; i++) H2Cmask[i] = (int *)malloc(H2Cmaskcolumn * sizeof(int));
    for (i = 0; i < RC1; i++) {
        for (j =0; j < N1; j++) {
            fscanf(fpr2, "%d", &H2Cmask[i][j]);
        }
    }
    fclose(fpr2);
   
    H1Crow = rc;
    H1Ccolumn = n;
    H1C = (int **)malloc(H1Crow * sizeof(int *));
    for (i = 0; i < H1Crow; i++) H1C[i] = (int *)malloc(H1Ccolumn * sizeof(int));
    n1 = n1 / 2;
    for (i = 0; i < H1Crow; i++) {
        for (j = 0; j < H1Ccolumn; j++) {
            H1C[i][j] = 0;
        }
    }
    int temp;
    for (i = 0; i < rc1; i++) {
        for (j = 0; j < n1; j++) {
            if (H1Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (H1Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e; m++) {
                   H1C[e * i + m][e * j + (m + H1Cmask[i][2 * j + 1]) % e] = 1; 
                }    
            }
               
        }
    }
    H2Crow = RC;
    H2Ccolumn = N;
    H2C = (int **)malloc(H2Crow * sizeof(int *));
    for (i = 0; i < H2Crow; i++) H2C[i] = (int *)malloc(H2Ccolumn * sizeof(int));
    N1 = N1 / 2;
    for (i = 0; i < H2Crow; i++) {
        for (j = 0; j < H2Ccolumn; j++) {
            H2C[i][j] = 0;
        }
    }
    for (i = 0; i < RC1; i++) {
        for (j = 0; j < N1; j++) {
            if (H2Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (H2Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e; m++) {
                   H2C[e * i + m][e * j + (m + H2Cmask[i][2 * j + 1]) % e] = 1; 
                }    
            }
               
        }
    }

    Hrow = e * rc1 + E * RC1;
    if (Hrow == (rc + RC)) printf("yes row = %d\n", Hrow);
    Hcolumn = e * n1 + E * N1;
    if (Hcolumn == (n + N)) printf("yes column = %d\n", Hcolumn);
    H = (int **)malloc(Hrow * sizeof(int *));
    for (i = 0;i < Hrow; i++) H[i] = (int *)malloc(Hcolumn * sizeof(int));

    H1row = Hrow;
    H1column = Hcolumn;
    H1 = (int **)malloc(H1row * sizeof(int *));
    for (i = 0; i < H1row; i++) H1[i] = (int *)malloc(H1column  * sizeof(int));

    for (i = 0; i < Hrow; i++) {
        for (j = 0; j < Hcolumn; j++) {
            H[i][j] = 0;
        }
    }
    //printf("123\n");
    printf("H1Crow = %d\n", H1Crow);
    for (i = 0; i < H1Crow; i++) {
        //printf("%d\n",i);
        //if (i = 1844) break;
        for (j = 0; j < H1Ccolumn; j++) {
            if (H1C[i][j] == 1){
                H[i][j] = H1C[i][j];
                H[rc + i][n + j] = H1C[i][j];
            }
        }
        //printf("j = %d\n", j);
        printf("a\n");
    }
    printf("123\n");
    for (i = 0; i < H2Crow; i++) {
        for (j = 0; j < H2Ccolumn; j++) {
            if (H2C[i][j] == 1) H[i][N + j] = H2C[i][j];
        }
    }
    int num = 0;
    int num1 = 0;
    int num2 = 0;
    printf("Hrow = %d; Hcolumn = %d\n", Hrow, Hcolumn);
    for (i = 0; i < Hrow; i++) {
        num = 0;
        for (j = 0; j < Hcolumn; j++) {
            printf("%d ", H[i][j]);
            H1[i][j] = H[i][j];
            if (H[i][j] == 1) num++;
        }
        printf("\n");
        printf("num = %d\n", num);
    }
    for (i = 0; i < Hcolumn; i++) {
        //num1 = 0;
        //num2 = 0;
        for (j = 0; j < Hrow; j++) {
            if (H[j][i] == 1) {
                if (i < n) num1++;
                else num2++;
            }
        }
        if (i < n) printf("num1 = %d", num1/(Hcolumn/2));
        else printf("num2 = %d", num2/(Hcolumn/2));
    }

    
    FILE *outfp1;
    outfp1 = fopen("C1parchematrix.txt","w");
    fprintf(outfp1,"%d ",Hcolumn);
    fprintf(outfp1,"%d ",Hrow);
    fprintf(outfp1,"\n");
    fprintf(outfp1,"%d ", num);
    fprintf(outfp1,"%d ", num1/(Hcolumn/2));
    fprintf(outfp1,"%d ", num2/(Hcolumn/2));
    fprintf(outfp1,"\n");
    for (i = 0; i < Hcolumn; i++) {
        for (j = 0; j < Hrow; j++) {
            if (H[j][i] == 1) fprintf(outfp1,"%d ",j+1);
        }
        fprintf(outfp1,"\n");
    }
    for (i = 0; i < Hrow; i++) {
        for (j = 0; j < Hcolumn; j++) {
            if (H[i][j] == 1) fprintf(outfp1,"%d ", j+1);
        }
        fprintf(outfp1, "\n");
    }
    fclose(outfp1);

    // GAUSS JORDAN METHOD
    int temprow = 0;    // for store pivot row 
    int stop = 0;
    int temptrans;
    int s = 0;
    printf("gauss jordan form!\n");
    for (i = 0; i < 3688 + s && temprow < 3688/*rc - s*//*n*/; i++) {                        // find pivot
        stop = 0;
        //printf("i = %d\n", i);
        for (j = temprow/*0*/; j < 3688 && stop == 0; j++) {      // where the pivot in the row
            if (H[j][i]  ==  1 && j == temprow) {
                temprow = j;
                //printf("j = %d temprow = %d \n", j, temprow);
                stop = 1;
            }
            else if (H[j][i]  ==  1 && j != temprow) {
                stop = 1;
                //printf("j = %d temprow = %d\n", j, temprow);
                for (k = 0; k < 7376; k++) {
                    temptrans = H[j][k];
                    H[j][k] = H[temprow][k];
                    H[temprow][k] = temptrans;
                }
            }
        }
        if(j == 3688 && stop == 0) {
            //i++;
            s++;
            printf("s = %d\n",s);
            continue;
        }
        // Eliminate the column have 1
        printf(" i = %d\n", i);
        for (j = 0; j < 3688; j++) {
            if (j == temprow) continue;                    // skip the original row because the row no need Elimination
            else {
                //printf("h\n");
                if (H[j][i] == 1) {
                    //printf("h\n");
                    for (k = 0; k < 7376; k++) {
                        H[j][k] = (H[j][k] + H[temprow][k]) % 2;
                    }
                }
            }
        }
        temprow++;
        //printf("temprow = %d\n", temprow);
    }
    int a = 0;
    int a1 = 0;
    int c[3688];
    int c1[3688];
    m = 0;
    int **Hsyst;
    int Hsystrow = 3688;
    int Hsystcolumn = 3688;
    int **Hsyst1;
    int Hsyst1row = 3688;
    int Hsyst1column = 3688;
    int **Gsys;
    int Gsysrow = 3688;
    int Gsyscolumn = 7376;
    int **G;
    int Grow = 3688;
    int Gcolumn = 7376;

    Hsyst = (int **)malloc(Hsystrow * sizeof(int *));
    for (i = 0; i < Hsystrow; i++) Hsyst[i] = (int *)malloc(Hsystcolumn * sizeof(int));
    Hsyst1 = (int **)malloc(Hsyst1row * sizeof(int *));
    for (i = 0; i < Hsyst1row; i++) Hsyst1[i] = (int *)malloc(Hsyst1column * sizeof(int));
    Gsys = (int **) malloc (Gsysrow * sizeof(int *));
    for (i = 0; i < Gsysrow; i++) Gsys[i] = (int *) malloc(Gsyscolumn *sizeof(int));
    G = (int **)malloc(Grow * sizeof(int *));
    for (i = 0; i < Grow; i++) G[i] = (int *)malloc(Gcolumn * sizeof(int));

    for (j = 0; j < 7376; j++) {
        m = 0;
        for (i = 0; i < 3688; i++) {
            if (H[i][j] == 0) m = m + 1;
        }
        if (m == 3688 - 1) {
            for (i = 0; i < 3688; i++) {
                Hsyst[i][a] = H[i][j]; 
            }
            c[a] = j;
            a++;
        }
        else {
            for (i = 0; i < 3688; i++) {
                Hsyst1[i][a1] = H[i][j];
            }
            c1[a1] = j;
            a1++;
        }
    }
    for (i = 0; i < 3688; i++) printf("%d ", c[i]);
    printf("\n");
    for (i = 0; i < 3688; i++) printf("%d ",c1[i]);
    printf("\n");
    m = 0;
    for (j = 0; j < 7376; j++) {
        //m = 0;
        if (j < 3688) {
            for (i = 0; i < 3688; i++) {
                Gsys[i][j] = Hsyst[i][j];
            }
        } else {
            /*for (i = 0; i < rc; i++) {*/
                for (k = 3688; k < 7376; k++) {
                    Gsys[m][k] = Hsyst1[k-3688][j-3688];
                }  
                m++;
            //}
            
        }
        
    }
    for (i = 0; i < 3688; i++) {
        for (j = 0; j < 7376; j++) {
            printf("%d ", Gsys[i][j]);
        }
    }
    //printf("nn\n");
    a = 0;
    a1 = 0;
    for (j = 0; j < 7376; j++) {
        if (c[a] == j) {
            printf("nna\n");
            for (i = 0; i < 3688; i++) {
                G[i][j] = Gsys[i][3688 + a];
            }
            a++;
        }
        if (c1[a1] == j) {
            printf("nna1\n");
            for (i = 0; i < 3688; i++) {
                G[i][j] = Gsys[i][a1];
            }
            a1++;
        }

    }
    printf("y\n");
    for (i = 0; i < 3688; i++) {
        for (j = 0; j < 7376; j++) {
            if (G[i][j] == 1) printf("%d ", j + 1);
        }
        printf("\n");
    }

    FILE *outfp;
    outfp = fopen("generator1.txt","w");
    for (int i = 0; i < 3688; i++) {
        for (int j = 0; j < 7376; j++) {
            fprintf(outfp,"%d ",G[i][j]);
        }
        fprintf(outfp,"\n");
    }
    fclose(outfp);

    FILE *outfp2;
    outfp2 = fopen("generator1.txt","w");
    for (int i = 0; i < 3688; i++) {
        printf("i = %d \n",i);
        for (int j = 0; j < 7376; j++) {
            if (G[i][j] == 1) fprintf(outfp2,"%d ",j+1);
            if (G[i][j] == 1) printf("%d ",j+1);
        }
        fprintf(outfp2,"\n");
    }
    fclose(outfp2);

    printf("\n\n");
    int temp1;
    for (i = 0; i < 3688; i++) {
        temp1 = 0;
        for (j = 0; j < 7376; j++) {
            temp1 += (G[i][j] * H1[i][j]);
            temp1 = temp1 % 2;   
        }
        printf("%d ", temp1);
        if(temp1 == 1) printf("yes\n");
    }


    for (i = 0; i < H1Cmaskrow; i++) free(H1Cmask[i]);
    free(H1Cmask);
    printf("a\n");
    for (i = 0; i < H2Cmaskrow; i++) free(H2Cmask[i]);
    free(H2Cmask);
    printf("a\n");
    for (i = 0; i < Hrow; i++) free(H[i]);
    free(H);
    printf("a\n");
    for (i = 0; i < H1row; i++) free(H1[i]);
    free(H1);
    printf("a\n");
    for (i = 0; i < H1Crow; i++) free(H1C[i]);
    printf("a\n");
    free(H1C);
    printf("a\n");
    for (i = 0; i < H2Crow; i++) free(H2C[i]);
    free(H2C);
    printf("a\n");
    for (i = 0; i < Hsystrow; i++) free(Hsyst[i]);
    free(Hsyst);
    for (i = 0; i < Hsyst1row; i++) free(Hsyst1[i]);
    free(Hsyst1);
    for (i = 0; i < Gsysrow; i++) free(Gsys[i]);
    free(Gsys);
    for (i = 0; i < Grow; i++) free(G[i]);
    free(G);

    return 0;

}