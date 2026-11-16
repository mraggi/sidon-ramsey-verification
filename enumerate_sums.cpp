// Independent enumeration via unordered PAIR SUMS, including diagonal pairs.
// No optimal-ruler length table is used. Cross-checks the difference enumerator.
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>
#include <omp.h>
using U=unsigned __int128;
struct B {U lo=0,hi=0;};
int n,target;bool normalized;
std::mutex mutex;
std::ofstream out;
std::atomic<uint64_t> count{0},nodes{0};
int pop(U x){return __builtin_popcountll(uint64_t(x))+__builtin_popcountll(uint64_t(x>>64));}
int first(U x){return uint64_t(x)?__builtin_ctzll(uint64_t(x)):64+__builtin_ctzll(uint64_t(x>>64));}
U right(B s,int x){return x?(s.lo>>x)|(s.hi<<(128-x)):s.lo;}
struct Search {
    int marks[128];uint64_t nn=0;
    std::vector<U> buffer;
    void flush(){std::lock_guard<std::mutex>lock(mutex);out.write(reinterpret_cast<char*>(buffer.data()),buffer.size()*16);count+=buffer.size();buffer.clear();}
    void dfs(int depth,U mask,B sums,U candidates){
        ++nn;
        if(depth==target){
            int translations=normalized?n-marks[depth-1]:1;
            for(int shift=0;shift<translations;++shift){buffer.push_back(mask<<shift);if(buffer.size()>=65536)flush();}
            return;
        }
        if(pop(candidates)<target-depth)return;
        while(candidates){
            int x=first(candidates);candidates&=candidates-1;
            B fresh{mask<<x,x?mask>>(128-x):0};
            if(2*x<128)fresh.lo|=U(1)<<(2*x);else fresh.hi|=U(1)<<(2*x-128);
            B next{sums.lo|fresh.lo,sums.hi|fresh.hi};
            U blocked=right(next,x);
            for(int i=0;i<depth;++i)blocked|=right(fresh,marks[i]);
            marks[depth]=x;
            dfs(depth+1,mask|(U(1)<<x),next,candidates&~blocked);
        }
    }
};
int main(int argc,char**argv){
    if(argc!=7){std::cerr<<"usage: sum_enumerate N SIZE THREADS full|BLOCKED_HI_HEX BLOCKED_LO_HEX OUT.bin\n";return 2;}
    n=std::stoi(argv[1]);target=std::stoi(argv[2]);int threads=std::stoi(argv[3]);
    if(n<1||n>127||target<1||target>14||threads<1)return 2;
    normalized=std::string(argv[4])=="full";
    U blocked=normalized?0:(U(std::stoull(argv[4],nullptr,16))<<64)|U(std::stoull(argv[5],nullptr,16));
    if(blocked>>n)return 2;
    U available=((U(1)<<n)-1)^blocked;
    out.open(argv[6],std::ios::binary);if(!out)return 2;
    auto start=std::chrono::steady_clock::now();
    if(normalized&&target>1){
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(int x=1;x<n;++x){
            Search s;s.marks[0]=0;s.marks[1]=x;
            B sums{1|(U(1)<<x),0};
            if(2*x<128)sums.lo|=U(1)<<(2*x);else sums.hi|=U(1)<<(2*x-128);
            U candidates=available&~((U(1)<<(x+1))-1);
            if(2*x<n)candidates&=~(U(1)<<(2*x));
            s.dfs(2,1|(U(1)<<x),sums,candidates);s.flush();nodes+=s.nn;
        }
    }else{
        // For one mark, direct enumeration avoids normalization duplicates.
        normalized=false;
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(int x=0;x<n;++x){
            if(!(available>>x&1))continue;
            Search s;s.marks[0]=x;B sums;
            if(2*x<128)sums.lo=U(1)<<(2*x);else sums.hi=U(1)<<(2*x-128);
            s.dfs(1,U(1)<<x,sums,available&~((U(1)<<(x+1))-1));s.flush();nodes+=s.nn;
        }
    }
    out.close();if(!out)return 1;
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::string meta="{\"n\":"+std::to_string(n)+",\"size\":"+std::to_string(target)+",\"count\":"+std::to_string(count.load())+",\"nodes\":"+std::to_string(nodes.load())+",\"seconds\":"+std::to_string(seconds)+",\"threads\":"+std::to_string(threads)+",\"complete\":true,\"method\":\"independent pair sums\"}";
    std::ofstream(std::string(argv[6])+".json")<<meta<<"\n";std::cout<<meta<<std::endl;
}
