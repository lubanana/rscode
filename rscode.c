#include<stdio.h>
#include <stdlib.h>  //define library
#include <math.h>
#include <string.h>

void galois(void);
int gf_mul(int x, int y);
int gf_add(int x, int y);
int inv_alpha(int x); 
int gf_division(int x, int y);
int * generator(void);
int * encoder(void);
int ** primitive(void);
bool codeword_check(void);

int prip[9][10] = { {0,}, {1,}, {2,}, 
                    {1,0,1,},        // 3
                    {1,0,0,1,},        // 4
                    {0,1,0,0,1,},    // 5
                    {1,0,0,0,0,1,},    // 6
                    {0,0,1,0,0,0,1,},    // 7
                    {0,1,1,1,0,0,0,1,} // 8
};
int prip3[3]={1,0,1}, prip4[4] = {1,0,0,1}, prip5[5] = {0,1,0,0,1}, prip6[6] = {1,0,0,0,0,1}; // , prip8[8]={0,1,1,1,0,0,0,1};
int m,t,n,k,powm;
int * alpha;
char * px;
int * gene_poly;
int * message;
int * en; 

void main(void)
{
    printf(" input the number n, k : ");
    scanf("%d,%d", &n, &k);
    frexp(n+1, &m);        // m = log2(n+1) 과 같은 결과
    m -= 1;                // 항상 값이 0이 되면 안되기 때문에 결과값에 +1이 된 상태이기 때문에 
                        // 제대로된 결과를 위해서 m + 1 
    powm = pow(2, m); // 2의 m승.
    t=(n-k)/2;
    galois();

    message = (int *)malloc(3 * sizeof(int));
    gene_poly = generator();
    message[0] = 5;
    message[1] = 3;
    message[2] = 1;
//    message[3] = 2;

    en = encoder();
    printf("The generator polynomial is : %d, %d, %d, %d, %d\n",gene_poly[0], gene_poly[1], gene_poly[2], gene_poly[3], gene_poly[4]);
    printf("The codeword polynomial is : %d   ,%d   ,%d   ,%d   ,%d   ,%d   ,%d   ",en[0],en[1],en[2],en[3],en[4],en[5],en[6]);
    bool result;
    result = codeword_check();
    printf("\nIs the codeword correct?  %d\n", result);  // 맞으면 1 출력.
    int ** irre = primitive();
    int x;
    int a = pow(2, m-1);
    for (x = 0; x < a; x++)
        printf("%d, %d, %d, %d, %d, %d\n", irre[x][0], irre[x][1], irre[x][2], irre[x][3], irre[x][4], irre[x][5]);
    



    char filename[7] = "mx.txt"; // m1.txt
    char temp[1];
    itoa(m, temp, 10);
    filename[1] = temp[0];
    FILE * f =fopen(filename,"w");
    int xc, yc;
    fprintf(f, "m is :%d\n", m); 
    fprintf(f, "덧셈\n");
    fprintf(f, "       ");
    for (xc = 0; xc < (powm - 1); xc++) {
        if (xc <= 9)
            fprintf(f, " a%d    ", xc);
        else
            fprintf(f, "a%d    ", xc);
    }
    fprintf(f, "\n");
    for (xc = 0; xc < (powm -1); xc++) {
        if (xc <= 9)
            fprintf(f, " a%d    ", xc);
        else
            fprintf(f, "a%d    ", xc);
        for (yc = 0; yc < (powm - 1); yc++) {
            if (gf_add(xc, yc) == -1)
                fprintf(f, "  %d    ", 0);
            else {
                if (gf_add(xc, yc) <= 9)
                    fprintf(f, " a%d    ", gf_add(xc, yc));
                else
                    fprintf(f, "a%d    ", gf_add(xc, yc));
            }
        }
        fprintf(f, "\n");
    }
    
    fprintf(f, "\n");
    fprintf(f, "곱셈\n");
    fprintf(f, "       ");
    for (xc = 0; xc < (powm - 1); xc++) {
        if (xc <= 9)
            fprintf(f, " a%d    ", xc);
        else
            fprintf(f, "a%d    ", xc);
    }
    fprintf(f, "\n");
    for (xc = 0; xc < (powm -1); xc++) {
        if (xc <= 9)
            fprintf(f, " a%d    ", xc);
        else
            fprintf(f, "a%d    ", xc);
        for (yc = 0; yc < (powm - 1); yc++) {
            if (gf_mul(xc, yc) == -1)
                fprintf(f, "  %d    ", 0);
            else {
                if (gf_mul(xc, yc) <= 9)
                    fprintf(f, " a%d    ", gf_mul(xc, yc));
                else
                    fprintf(f, "a%d    ", gf_mul(xc, yc));
            }
        }
        fprintf(f, "\n");
    }


    /*int x;
    int y;
    for (x = -1; x < 7; x++) {
        for (y = -1; y < 7; y++) {
            printf("x : %d y : %d 결과 : %d\n", x, y, gf_add(x, y));
        }
        printf("\n");
    }    */    
}

/*alpha table. */
void galois(void) 
{    
    alpha = (int *)malloc(powm * sizeof(int)); //powm 만큼 배열 만들기.
    px = (char *)malloc((m + 1) * sizeof(char));// 배열 선언시, 변수 사용불가.그러므로 malloc 사용.
    int i,x,y,pinit;
    
    pinit = 0;
    memset(px, 0, (m + 1) * sizeof(char)); // px 초기화
    px[1] = 1;  //  alpha 0 의 상태로 초기화 즉, 첫번째 자리에 1 넣고, 나머지는 px초기화에 의해 0이 들어가 있음.
    alpha[0] = 1; // alpha0 =1
    alpha[n] = alpha[0];  // finit field의 성질.
    printf("%5s%5s%5s%5s%5s%5s%5s%5s%5s%10s%10s\n", "i","p1","p2","p3","p4","p5","p6","p7","p8","pinit","alpha[i]");
     // for문 밖에 있어야 함.
    for (i = 1; i < n; i++)   //alpha[1] ~ alpha[n-1]  
    {    
        pinit = px[m];
        for (x = m - 2; x >= 0; x--) 
        {
            if (prip[m][x] == 0)
                px[x + 2] = px[x + 1];
            else 
                px[x + 2] = px[x + 1] ^ pinit;
        }
        px[1] = pinit; // 여기까지가 각 alpha[]내부에 수를 넣는 일. 
        
        alpha[i] = 0; //  반드시, 해야함.초기화.
        for (y = 1; y <= m; y++)  
            alpha[i] += px[y] * (int)pow(2, y-1); //alpha[i]에 십진수로 저장. 메모리 상에는 2진수로 저장됨.
    
        printf("%5d%5d%5d%5d%5d%5d%5d%5d%5d%10d%10d\n",i,px[1],px[2],px[3],px[4],px[5],px[6],px[7],px[8],pinit,alpha[i]);
    }
}

/* multiplication in the Galois Field*/
int gf_mul(int x, int y) // alpha의 지수값을 받기. return 값도 지수값으로. 
{
    if(x== -1 || y== -1) return 0; // 나중에 지수값을 대입시, 숫자 0은 -1 대입해야 한다. 
    else return (x + y) % n;
}

/*addition in the Galois Field*/
int gf_add(int x, int y) // 지수 값을 받아서 지수값을 돌려준다. 
{
    int xtemp, ytemp, temp, z;
    xtemp = alpha[x];
    ytemp = alpha[y];
    if (x == -1)
        xtemp = 0;
    if (y == -1)
        ytemp = 0;
    temp = xtemp ^ ytemp;    // x와 y의 값을 xor한 결과를 temp 에 저장 한 뒤에
    for (z = 0; z < n; z++) 
    {
        if (temp == alpha[z])    // 테이블에서 temp와 같은 값을 찾아서 
            return z;            // 그 값에 해당되는 지수를 리턴한다. 
    }
    return -1; // 자기자신끼리 합 = 0  0+0 =0 (=-1 in this progrma), alpha[i]중에서 그 값을 찾을수가 
                // 없으므로, return 값을 -1로 하였다. 
}

int inv_alpha(int x) // 역수 구하는 것.
{
    return n-x; // GF에서 서로 곱해 1이 되는 것은 두 지수의 합이 n일때 이므로.
}

/*division in the GF*/
int gf_division(int x, int y) 
{
    return gf_mul(x, inv_alpha(y));    // 나눗셈은 x와 y의 역수가 곱해진 결과 
}



/* generating polynomial mutiplication*/
int * generator(void) 
{
    int * gene_poly = (int *)malloc((2*t + 1) * sizeof(int));  //일종의 임시배열이지만, 나중에 최종결과값 저장할 배열.
    int * temp = (int *)malloc((2*t + 1) * sizeof(int)); //계산상 필요한 임시 배열.
    int gp[2]; // 일차다항식의 계수를 받기 위한 배열. 
    int i,x,y,z, maxdegree;
    
    memset(gene_poly, -1, (2*t+1) * sizeof(int)); // gene_poly 배열을 -1 (=숫자 0) 으로 초기화.
    
    gp[1] = 0; // X의 계수인 1을 alpha[1]로 읽지 않게 하기 위해 0, 즉 alpha[0]으로 초기화.
    /* [X+alpha] 를 temp에 대입하기. */
    temp[0] = 1; 
    temp[1] = gp[1];
    /*한 항씩 곱해가면서 temp에 값을 넣음 */
    for (x=    0; x < 2*t-1; x++)     // 2*t-1번의 계산 필요. 
    {    
        maxdegree = x + 2;  //temp에 들어 있는 최고차항의 수.  (한번 계산할때 2차가 됨을 떠올리면 됨)
        for (y=0; y < 2; y++)   //일차다항식만 계속해서 곱해지므로, 두번만 돌리면 됨. 
        {
            gp[0] = x + 2;   // 상수항에는 alpha의 제곱, 세제곱... 등의 값이 들어감을 표현, 알파의 지수승만
                                 //따야 하므로 x+2와 같이 표현.
            for (z=0; z < maxdegree; z++) // temp에 들어있는 항의 수만큼 z이 바뀜.
            {
                if (gene_poly[y + z] != -1) // gene_poly에 어떠한 값이 미리 들어가 있으면, 앞서, 그 차수에서 계산 되었단
                                               //뜻이므로 거기에다가 더하고,
                    gene_poly[y + z] = gf_add(gene_poly[y + z], gf_mul(gp[y], temp[z]));    
                else  //  gene_poly에 어떤 값도 없으면,  그냥 곱한 값을 넣는다.  
                    gene_poly[y + z] = gf_mul(gp[y], temp[z]);
            }
            
        }
        if (x < 2 * t - 2)    // 마지막 계산이 아닐 때는 gene_poly 초기화, temp 에 다음계산할 식 대입
        {
            for (i = 0; i < (2 * t + 1); i++)
            {
                    temp[i] = gene_poly[i];  //gene_poly 에 있는 값을 temp에 집어 넣는다. 
                    gene_poly[i] = -1;  // gene_poly 를 비운다.   0으로 초기화 시킴. 
            }
        }
        else    // 계산이 끝났을 때는 gene_poly 초기화와 temp에 다음 계산 식 대입이 필요없으니 
            return gene_poly;    // gene_poly 에 들어있는 계산 결과를 리턴해준다 
    }
}

int * encoder(void) 
{
    int * en = (int *)malloc(n * sizeof(int));  // encoding 한 결과를 넣을 동적배열. 
    int * regi = (int *)malloc((2 * t) * sizeof(int));  // register 를 2t 개 만큼 만듬
    int x, y, i, feedback_temp, feedback;

    // 가장 먼저 들어오는 메시지는 0 번
    // 인코딩 결과 배열에 미리 메시지 쉬프트 시켜서 입력.
    for (x = 0; x < k; x++) 
        en[n - (x + 1)] = message[x];

    memset(regi, -1, (2 * t) * sizeof(int)); // regi 초기화 
    
    feedback = -1;  // feedback 초기화
    for (y = 0; y < k; y++)  // 전체 clock 수 = k
    {  
        if (feedback == -1) // 첫번째 clock 을 의미함.
        {   // 메시지 *  gene_poly 계수             
            for (x = 0; x < (2 * t); x++) 
                regi[x] = gf_mul(message[y], gene_poly[x]); // 메시지 심벌도 지수만 받도록 짜야함
        }
        else // 1클럭이 지난 후에는 여기로 
        {    
            for (x = 0; x < (2 * t); x++)
            {
                if(regi[x] != -1)  // regi에 수가 들어있으면 그 수에다가 [feedback * gene_poly 각 계수]를 더함. 
                    regi[x] = gf_add(regi[x],gf_mul(feedback, gene_poly[x]));
                else  // regi에 수 없으면, 그냥 곱한것을 넣음.
                    regi[x] = gf_mul(feedback, gene_poly[x]);
            }
        }
        if (y < k-1)  // 쉬프트는 k-1 번 만큼 시키면 되므로, (=마지막에서는 쉬프트 시킬 필요 없음)
        {
            feedback_temp = regi[2 * t - 1]; //register 가장 마지막 값을 임의로 넣어둠 
            // 쉬프트 
            for (x = (2 * t) - 1; x >= 0; x--) 
                regi[x] = regi[x - 1];
            regi[0] = -1;  // 숫자 0을 넣는 것임 , 즉 비워두기(초기화). 
            feedback = gf_add(feedback_temp, message[y + 1]); // feedback 값은 다음 message와 그전에 저장시킨 feedback_temp를 더해서 구함. 
        }
    }
    for (i=0; i < (2*t); i++) // temp에 저장된값을 en에 넣기. 
        en[i] = regi[i];
    
    return en;
}

/* codeword check*/
bool codeword_check(void)
{
    int x,y,temp;
    temp = -1; // temp 초기화.

    for (x = 0; x < 2*t; x++)
    {
        for(y = 0; y < n; y++)
        {
             gf_add(temp, gf_mul(en[y],gene_poly[x]));    // temp에 en[y] * gene_poly[x] 대입
        }
    }

    if (temp != -1)    // temp에 들어있는 값이 0이 아니라면 인코딩 된 식이 잘못된 것
        return false;
    return true;    // 들어있는 값이 0이면 인코딩된 식이 제대로 된 것 
}

/*primitive polynomial*/
int ** primitive(void)
{
    int x, y, z, bin, inv, temp, counter, feedback, feedback_temp;
    int a = pow(2, m-1);
    int ** irre_poly = (int **)malloc( a * sizeof(int *));  //a 층 만큼 만들어놓기 
    /* condition1 for the primitive polynomial*/
    for (x = 0; x < a; x++)
    { 
        irre_poly[x] = (int *)malloc( (m+1) *sizeof(int));  // 각 층에 들어가는 배열 수.
        irre_poly[x][0] = irre_poly[x][m] = 1; // 배열의 첫번째 자리와 마지막 자리에는 항상 1을 대입.
        bin = x;
        inv = m - 2;
        temp = 0;
        for (y = 1; y < m; y++)
        {
            if (bin & (int)pow(2, inv)) // 비트 연산시 자료형을 맞춰주기 위해 bin의 int형과 같이pow()도 int형으로 해줌.
                irre_poly[x][y] = 1; // if 괄호 안이 true이면(=수가 들어가 있으면), irre_poly[x][y]=1 수행.  
            else
                irre_poly[x][y] = 0;
            inv--;
            temp += irre_poly[x][y]; // 첫자리와 마지막 자리 제외하고 가운데 수들을 더함. 
        }
        if(temp%2 == 0) // 그 더한수가 짝수일 경우,  편의상 -1을 적어둬서 나중에 조건2를 행할때 제외시킴.
            irre_poly[x][0] = -1;
    }    
    /* condition 2 for primitive polynomial */
    int * remainder = (int *)malloc( m * sizeof(int)); 
    memset(remainder, 0, m * sizeof(int));
    feedback = -1;
    for (bin = 0; bin < a; bin++) 
    {
        if (irre_poly[bin][0] == -1)
            continue;
        for (x = m +1; x < n ; x++)
        {
            inv = x + 1;
            for (counter = 0; counter < inv; counter++)  // temp 위해서. X^4+1,X^5+1..발생.
            {  
                if (counter == 0 || counter == inv - 1)
                    temp = 1;
                else
                    temp = 0;
            
                if(feedback ==-1)      // 첫번째 clock 의미.
                { 
                    for(z =0; z < m; z++) 
                        remainder[z] = irre_poly[bin][z] *    temp;
                }
                else         // 첫번째 clock 이후. 
                { 
                    for(z =0; z < m; z++)
                    {
                        if (remainder[z] != 0)
                            remainder[z] ^= feedback * irre_poly[bin][z]; 
                        else 
                            remainder[z] = feedback * irre_poly[bin][z];
                    }
                }

                if (counter < m + 1)
                {
                    feedback_temp = remainder[m - 1];
                    // 쉬프트
                    for (z = m - 1; z; z--)
                        remainder[z] = remainder[z - 1];
                    remainder[0] = 0;
                    if (counter + 1 != m + 1)
                        temp = 0;
                    else
                        temp = 1;
                    feedback = feedback_temp ^ temp;
                }
            }
            
            temp = 0;
            for (z = 0; z < m; z++) 
                temp += remainder[z];
            if (temp == 0) 
            {
                irre_poly[bin][0] = -1;
                x = n;
                //continue;
            }
            for (z = 0; z < m; z++)
                remainder[z] = 0;
            feedback = -1;
            feedback_temp = 0;
        }
        
    }

    return irre_poly;
}
