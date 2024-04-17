#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"
using namespace std;
using namespace emp;


const static int nP = 8;
const static int BENCH = 10;
int party, port;

// const char out3[] = "7e50ebb4bd718d7b";

int main(int argc, char** argv) {
 
	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
    if(party > nP)return 0;

    struct timeval t1,t2;
    double timeuse = 0;

    // string file = "/home/zck/文档/MP-ABY/emp-agmpc/test/adder64.txt";
    string file = "/home/zck/文档/Esorics-MPABY/MPABY/test/adder64.txt";
	BristolFormat cf(file.c_str());
    // BristolFormat cf(file2.c_str());
	

    NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(nP);	
    ThreadPool pool2(nP);	
    PRG prg;

    CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);

    MPGC<nP>* GC = new MPGC<nP>(ios, &pool2, party, &cf);

    ShareA r, res, r_A;
    ShareB a;
    bool a_value,c;
    int64_t rec_value = 0, r_value = 0;
  
    ShareB res_vec[65];
     ShareB res_vec1[65];
    prg.random_bool(&a_value,1);
    prg.random_data(&r_value,sizeof(int64_t));

    mpc->GMW_A->set_ShareA(&r,r_value,2);


    A2B_ctx2 ctx;
    A2B_ctx2 ctx1;

    mpc->A2B_setup_api(&ctx,&r,GC);



    gettimeofday(&t1,NULL);
    mpc->A2B_online_api(res_vec,&ctx,&r,GC);

    gettimeofday(&t2,NULL);
    timeuse =  timeuse + (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        

    bool rec_revc[64];
    mpc->GMW_B->rec_ShareB(rec_revc,res_vec,64);


    string res_string = "";
    for (int i = 63; i >= 0; i--)
    {
        res_string += (rec_revc[i]?"1":"0");
        cout<<rec_revc[i];

    }
    cout<<endl;
        
    mpc->GMW_A->rec_ShareA(rec_value,&r);

    char out3[100] = {0};
    sprintf(out3,"%16lx",rec_value);
    cout<<hex_to_binary(string(out3))<<endl;

    cout <<"1:"<< (res_string == hex_to_binary(string(out3))? "GOOD!":"BAD!")<<endl<<flush;

    delete mpc;
    // cout<<"A2B_online throughput: "<<BENCH/timeuse<<"opt/s"<<endl;
    cout<<"A2B_online runtime: "<<(timeuse/1)*1000.0<<"ms"<<endl;


}