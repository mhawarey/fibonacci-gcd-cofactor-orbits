// paper7_search_original.c: the program that wrote search_logs/ (Theorem 6). Starts: u = alpha+1, v = beta+1 with v - u*phi in (-1/phi, 1/phi^2).
// Finds every cycle of R whose minimal-alpha state has alpha in [A0,A1], all states < M (as __int128), period<=CAP.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128; typedef unsigned long long u64;
static u128 gcd128(u128 x,u128 y){while(y){u128 t=x%y;x=y;y=t;}return x;}
static u128 isqrt128(u128 n){ u128 x=(u128)__builtin_sqrtl((long double)n); while(x*x>n)x--; while((x+1)*(x+1)<=n)x++; return x;}
// floor(u*phi) = floor((u + sqrt(5u^2))/2)
static void pr(u128 x){char b[50];int i=49;b[i]=0;if(!x){printf("0");return;}while(x){b[--i]='0'+(int)(x%10);x/=10;}printf("%s",b+i);}
int main(int argc,char**argv){
  u64 A0=strtoull(argv[1],0,10),A1=strtoull(argv[2],0,10); int CAP=atoi(argv[3]);
  u128 M=((u128)1)<<100; u64 tested=0;
  for(u64 al=A0;al<=A1;al++){
    u128 u=al+1; u128 s=isqrt128(5*u*u); u128 v=(u+s)/2; // floor(u phi) (s=floor(u sqrt5), exact since sqrt5 irrational)
    // unique v with v - u*phi in (-1/phi, 1/phi^2): v=floor(u phi) if frac<1/phi else floor+1
    { u128 x=u-1, y=v-1; if(!(5*x*x < (2*y-x)*(2*y-x))) v=v+1; }
    u128 a=al, b=v-1; if(b<1) continue; if(gcd128(a,b)!=1) continue; tested++;
    u128 X=a,Y=b; int k=0,i;
    for(i=1;i<=CAP;i++){ u128 g=gcd128(Y,X+1); u128 nx=Y/g, ny=(X+Y+1)/g; X=nx;Y=ny; if(g>1)k++;
      if(Y>M){i=-1;break;}
      if(X<a){i=-1;break;}
      if(X==a&&Y==b)break; }
    if(i>0&&i<=CAP){ printf("%d %d ",i,k); pr(a); printf(" "); pr(b); printf("\n"); fflush(stdout);}
  }
  fprintf(stderr,"tested %llu\n",tested); return 0;}
