// paper7_search.c: the search of Theorem 6 (Hawarey 2026, "A Size Law for Periodic Orbits of the Fibonacci-gcd Cofactor Map").
// For alpha in [A0, A1]: the golden-line start of Lemma 1 (and, as a margin, the other candidate); iterate R in 128-bit
// integers; print "alpha beta steps amplifications reason" with reason 0 = returned (a cycle), 1 = alpha dropped,
// 2 = an entry above 2^100, 3 = 5,000 steps. Usage: ./paper7_search A0 A1 5000 [--cycles-only]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned __int128 u128; typedef unsigned long long u64;
static u128 gcd128(u128 x,u128 y){while(y){u128 t=x%y;x=y;y=t;}return x;}
static u128 isqrt128(u128 n){ u128 x=(u128)sqrtl((long double)n); while(x*x>n)x--; while((x+1)*(x+1)<=n)x++; return x;}
int main(int argc,char**argv){
  u64 A0=strtoull(argv[1],0,10),A1=strtoull(argv[2],0,10); int CAP=atoi(argv[3]); u128 M=((u128)1)<<100; int only=(argc>4); u64 tested=0;
  for(u64 al=A0;al<=A1;al++){
    u128 u=al+1; u128 s=isqrt128(5*u*u); u128 v=(u+s)/2;
    { u128 x=u-1, y=v-1; if(!(5*x*x < (2*y-x)*(2*y-x))) v=v+1; }
    u128 a=al, b=v-1; if(b<1) continue; if(gcd128(a,b)!=1) continue;
    u128 X=a,Y=b; int k=0,i,reason=3;
    for(i=1;i<=CAP;i++){ u128 g=gcd128(Y,X+1); u128 nx=Y/g, ny=(X+Y+1)/g; X=nx;Y=ny; if(g>1)k++;
      if(Y>M){reason=2;break;} if(X<a){reason=1;break;} if(X==a&&Y==b){reason=0;break;} }
    if(!only || reason==0) printf("%llu %llu %d %d %d\n",(u64)a,(u64)b,i>CAP?CAP:i,k,reason);
    tested++;
  }
  fprintf(stderr,"tested %llu\n",tested); return 0;}
