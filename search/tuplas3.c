// tuplas N SS11.txt R OUT : list disjoint R-tuples (i<j<..) of 11-sets, canonical up to x->N+1-x
// (keep the tuple if its sorted index tuple <= the sorted tuple of reflected indices).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
int n, N, R, P=1, PR=0; u128 *M; int *ref; FILE *out; long tot, can;
int cur[8];
static int cmpi(const void*a,const void*b){return *(int*)a-*(int*)b;}
void rec(int k, int from, u128 used) {
    if (k == R) {
        tot++;
        int r[8]; for (int i=0;i<R;i++) r[i]=ref[cur[i]];
        qsort(r,R,sizeof(int),cmpi);
        for (int i=0;i<R;i++) { if (cur[i]<r[i]) break; if (cur[i]>r[i]) return; }
        can++;
        if (out) { for (int i=0;i<R;i++) fprintf(out,"%d%c",cur[i],i<R-1?'_':'\n'); }
        return;
    }
    for (int j=from;j<n;j++) if (!(M[j]&used)) { if (k==0 && j%P!=PR) continue; cur[k]=j; rec(k+1,j+1,used|M[j]); }
}
int main(int argc,char**argv){
    N=atoi(argv[1]); R=atoi(argv[3]); FILE*f=fopen(argv[2],"r"); char line[256];
    M=malloc(60000*sizeof*M); ref=malloc(60000*sizeof(int)); u128 *RM=malloc(60000*sizeof*RM);
    while(fgets(line,sizeof line,f)){ u128 m=0,rm=0; char*p=line; while(*p){ long x=strtol(p,&p,10); if(x){m|=(u128)1<<x; rm|=(u128)1<<(N+1-x);} if(*p==',')p++; else if(*p=='\n')break; } M[n]=m; RM[n]=rm; n++; }
    for(int i=0;i<n;i++){ ref[i]=-1; for(int j=0;j<n;j++) if(M[j]==RM[i]){ref[i]=j;break;} if(ref[i]<0){fprintf(stderr,"no reflection for %d\n",i);return 1;} }
    out = argc>4 ? (strcmp(argv[4],"-")?fopen(argv[4],"w"):stdout) : NULL; if(argc>6){P=atoi(argv[5]);PR=atoi(argv[6]);}
    rec(0,0,0);
    fprintf(stderr,"N=%d R=%d sets=%d tuples=%ld canonical=%ld\n",N,R,n,tot,can);
}
