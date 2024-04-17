#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"

using namespace std;
using namespace emp;


const static int nP = 3;
int party, port;

int main(int argc, char** argv) {
 
	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
	if(party > nP)return 0;

	NetIOMP<nP> io(party, port);
	ThreadPool pool(4);	
    PRG prg;

    NetIOMP<nP> io_1(party, port+2*(nP+1)*(nP+1)+2);
	NetIOMP<nP> io_2(party, port+2*(nP+1)*(nP+1)+3);
	NetIOMP<nP> *ios[2] = {&io_1, &io_2};
	// string file = circuit_file_location+"/AES-non-expanded.txt";
	// file = circuit_file_location+"/sha-1.txt";
	string file = "/home/jiabin/extdisk/zck/Documents/MPC/MP-ABY/emp-agmpc/test/adder64.txt";
	BristolFormat cf(file.c_str());
	


    CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);
    cout <<"Setup MPC:\t"<<party<<"\n";

    ShareB in[128];
    ShareB out[128];
    bool in1_value[128] = {0};

    bool in_delta[128] = {0};
    bool in_Public[128] = {0};

    bool out_delta[128] = {0};
    bool out_delta1[128] = {0};
    bool out_Public[128] = {0};

    in1_value[0] = 1;
    in1_value[64] = 1;

    for (int k = 0; k < 100; k++)
    {

    

    for (int i = 0; i < 128; i++)
    {
           mpc->GMW_B->set_ShareB(in+i,in1_value[i],1);
           in_delta[i] = in[i].delta;
           in_Public[i] = in[i].PublicVal;
           bool rec;
           mpc->GMW_B->rec_ShareB(rec,in+i);
        //    cout<<"rec:"<<rec<<"\tvalue:"<<in1_value[i]<<endl;
    }

    // CMPC<nP>* GC = new CMPC<nP>(ios, &pool, party, &cf,in_delta);
    MPGC<nP>* GC = new MPGC<nP>(ios, &pool, party, &cf);

	cout <<"Setup GC:\t"<<party<<"\n";

    GC->function_independent(in_delta);
	cout <<"FUNC_IND:\t"<<party<<"\n";	

    

    GC->function_dependent(out_delta);
	cout <<"FUNC_DEP:\t"<<party<<"\n";

    GC->online(in_Public, out_Public);

    // for (int i = 0; i < cf.n3; i++)
    // {
    //     cout << "out_delta:" << out_delta1[i] << " out_delta:" << out_delta[i]<< endl;
    //     /* code */
    // }
    

    for (int i = 0; i < cf.n3; i++)
    {
        out[i].delta = out_delta[i];
        out[i].PublicVal = out_Public[i];
        bool res;
        mpc->GMW_B->rec_ShareB(res,out+i);
        cout<<res;
    }
    cout<<endl;

    delete GC;
            /* code */
    }
    

	// uint64_t band2 = io_1.count();
	// cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
	// cout <<"ONLINE:\t"<<party<<"\n";

    
}