#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"
// #include "../../MPABY_GC/GC_mpc.h"
using namespace std;
using namespace emp;

#define test_B2A
#define BENCH 1
const static int nP = 32;
int party, port;

// const char out3[] = "7e50ebb4bd718d7b";

int main(int argc, char** argv) {
 
	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
	if(party > nP)return 0;

	NetIOMP<nP> io(party, port);
	ThreadPool pool(nP);	
    PRG prg;

    NetIOMP<nP> io_1(party, port+2*(nP+1)*(nP+1)+2);
	NetIOMP<nP> io_2(party, port+2*(nP+1)*(nP+1)+3);
	NetIOMP<nP> *ios[2] = {&io_1, &io_2};


    struct timeval t1,t2;
    double timeuse = 0;
    

#ifdef test_Bit2A

    for (int i = 0; i < 10; i++)
    {  
        CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);
        // cout <<"Setup MPC:\t"<<party<<"\n";

        // cout <<"Setup GC:\t"<<party<<"\n";  
        ShareA res, r_A;
        ShareB a;
        bool a_value;
        int64_t rec_value = 0;

        bool d_delta[64] = {0};
        ShareB ri[64];
        ShareB *res_vec = new ShareB[65];
        prg.random_bool(&a_value,1);
        // prg.random_data(&r_value,sizeof(int64_t));
        // r_value = 0x7e50ebb4bd718d7b;
        mpc->GMW_B->set_ShareB(&a,a_value,3);
        // mpc->GMW_A->set_ShareA(&r,r_value,1);

        // mpc->GMW_A->mulA_constant(&res,&r,8589934592);
        Bit2A_ctx2 ctx;
        mpc->Bit2A_setup_api(&ctx,a.delta);
        mpc->Bit2A_online_api(&res,&ctx,&a);
        // mpc->Bit2A_setup(c,&r,&a);
        // mpc->Bit2A_online(&res,&a,c,&r);
        mpc->GMW_A->rec_ShareA(rec_value,&res);
        cout<<"delta:"<<res.delta<<"\tDelta:"<<res.PublicVal<<endl;
        printf("res:%d ",rec_value);
        printf("res_a:%d\n",a_value);
        cout <<((rec_value != (int64_t)a_value)? "fail":"success")<<endl;

        delete mpc;


    }

#endif
    
#ifdef test_B2A

        CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);
        // cout <<"Setup MPC:\t"<<party<<"\n";

        ShareA res, r_A;
        ShareB a[64];
        bool a_value[64];
        bool a_rec[64] = {0};
        int64_t rec_value = 0;
        int64_t d = 0;
        bool a_delta[64] = {0};

        prg.random_bool(a_value,64);
        for (int i = 0; i < 64; i++)
        {
            mpc->GMW_B->set_ShareB(a+i,a_value[i],3);
            a_delta[i] = a[i].delta;
            mpc->GMW_B->rec_ShareB(a_rec[i],a+i);
        }
    
        Bit2A_ctx2 ctx[64];
        int64_t res_delta = 0;
        

        mpc->B2A_setup_api(ctx,res_delta, a_delta);




        gettimeofday(&t1,NULL);
        mpc->B2A_online_api(&res,ctx,res_delta,a);

        gettimeofday(&t2,NULL);
        timeuse =  timeuse + (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        

        mpc->GMW_A->rec_ShareA(rec_value,&res);

        if (party == 1)
        {
            string res_string = "";
            for (int i = 63; i >= 0; i--)
            {
                res_string += (a_rec[i]?"1":"0");
                cout<<a_rec[i];
            }
            cout<<endl;
            char out3[100] = {0};
            sprintf(out3,"%16lx",rec_value);
            cout<<hex_to_binary(string(out3))<<endl;

            cout <<"1:"<< (res_string == hex_to_binary(string(out3))? "GOOD!":"BAD!")<<endl<<flush;
            // cout<<"B2A_online throughput: "<<BENCH/timeuse<<"opt/s"<<endl;
            cout<<"B2A_online runtime: "<<(timeuse/BENCH)*1000.0<<"ms"<<endl;
        }
        delete mpc;
        

#endif

}