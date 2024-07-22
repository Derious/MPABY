#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"
using namespace std;
using namespace emp;

#define test_ReLU
#define _debug
const static int nP = 16;
const static int BENCH = 1;
const static int length = 1;
int party, port;

// const char out3[] = "7e50ebb4bd718d7b";
int64_t ReLU(int64_t x) {
    return x > 0 ? x : 0;
}

int main(int argc, char** argv) {
 
	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
    if(party > nP)return 0;

    struct timeval t1,t2;
    double timeuse = 0;

    struct timeval total1,total2;
    double timeuse_total = 0;

    string file = "/home/zck/文档/MPABY/MPABY/test/adder64-fushion.txt";
	BristolFormat cf(file.c_str());
    // BristolFormat cf(file2.c_str());
	

    NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(nP);	
    ThreadPool pool2(nP);	
    PRG prg;

    CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);

    #ifdef test_A2B
    vector<MPGC<nP>*> GC;

    for (int i = 0; i < length; i++)
    {
        GC.push_back(new MPGC<nP>(ios, &pool2, party, &cf));
    }
    
    ShareA res, r_A;
    ShareA* r = new ShareA[length];
    // ShareB a;
    bool a_value;
    
    int64_t* rec_value = new int64_t[length];
  
    ShareB res_vec[length][65];
    // ShareB res_vec1[65];
    prg.random_bool(&a_value,1);
    
    for (uint16_t i = 0; i < length; i++)
    {
        int64_t r_value = 0;
        prg.random_data(&r_value,sizeof(int64_t));
        mpc->GMW_A->set_ShareA(&r[i],r_value,2);
    }
    
    


    A2B_ctx* ctx;
    ctx = new A2B_ctx[length];

    for (uint16_t i = 0; i < length; i++)
    {
        mpc->A2B_setup_api(&ctx[i],r[i].delta,GC[i]);
    }


    for (uint16_t i = 0; i < length; i++)
    {
        mpc->A2B_online_api(res_vec[i],&ctx[i],&r[i],GC[i]);
    }

    bool rec_revc[length][64];
    for (uint16_t i = 0; i < length; i++)
    {
        mpc->GMW_B->rec_ShareB(rec_revc[i],res_vec[i],64);
    }
    // mpc->GMW_B->rec_ShareB(rec_revc,res_vec,64);


    string res_string[length];
    string res_string1 = "";
    for (uint16_t j = 0; j < length; j++)
    {
       

        mpc->GMW_A->rec_ShareA(rec_value[j],&r[j]);

        if(party == 1) {
            for (int i = 63; i >= 0; i--)
            {
                res_string[j] += (rec_revc[j][i]?"1":"0");
                cout<<rec_revc[j][i];

            }
            cout<<endl;
            char out3[100] = {0};
            sprintf(out3,"%16lx",rec_value[j]);
            cout<<hex_to_binary(string(out3))<<endl;
            cout <<"1:"<< (res_string[j] == hex_to_binary(string(out3))? "GOOD!":"BAD!")<<endl<<flush;
        }
    }

  
    gettimeofday(&t1,NULL);
    for (int i = 0; i < BENCH; i++)
    {
        mpc->A2B_online_api(res_vec[0],&ctx[0],&r[0],GC[0]); 
    }
    gettimeofday(&t2,NULL);
    timeuse =  timeuse + (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
    
        

    delete mpc;
    // cout<<"A2B_online throughput: "<<BENCH/timeuse<<"opt/s"<<endl;
    if (party == 1)
    {
       cout<<"A2B_online runtime: "<<(timeuse/BENCH)*1000.0<<"ms"<<endl;
    }

    #endif
    

    #ifdef test_ReLU
    for (int i = 0; i < BENCH; i++)
    {

        MPABY_ReLU_CTX ctx;
        MPGC<nP>* GC = new MPGC<nP>(ios, &pool2, party, &cf);


        int64_t a=0;
        ShareA a_maskA;
        prg.random_data(&a,sizeof(int64_t));
        a = a >> 1;


        mpc->GMW_A->set_ShareA(&a_maskA,a,1);

        ShareA res;

        gettimeofday(&total1,NULL);
        mpc->MPABY_ReLU_setup(ctx.delta_u,&(ctx.ctx),GC,ctx.t,&(ctx.r),ctx.delta_ra,a_maskA.delta);
        res.delta = ctx.delta_u;

    
        gettimeofday(&t1,NULL);
        mpc->MPABY_ReLU_online(res.PublicVal,ctx.delta_u,&(ctx.ctx),GC,ctx.t,&(ctx.r),ctx.delta_ra,a_maskA.PublicVal, a_maskA.delta);
        gettimeofday(&t2,NULL);
        gettimeofday(&total2,NULL);
        timeuse =  timeuse + (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        timeuse_total =  timeuse_total + (total2.tv_sec - total1.tv_sec) + (double)(total2.tv_usec - total1.tv_usec)/1000000.0;
        
        #ifdef _debug
        int64_t rec_A;
        int64_t res_v;
        mpc->GMW_A->rec_ShareA(rec_A,&a_maskA);
        mpc->GMW_A->rec_ShareA(res_v,&res);

        printf("a:%ld\n",rec_A);
        cout <<"res:"<<res_v<<" res:"<< ReLU(rec_A) << ((res_v == ReLU(rec_A))? " GOOD!":" BAD!")<<endl<<flush;
        #endif
        delete GC;
    }
    cout<<"MPABY_ReLU_online runtime: "<<(timeuse/BENCH)*1000.0<<"ms"<<endl;
    cout<<"MPABY_ReLU_total runtime: "<<(timeuse_total/BENCH)*1000.0<<"ms"<<endl;

    #endif


}