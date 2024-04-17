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

    int64_t r_value = 0;
    int64_t mask_input = 0;
    prg.random_data(&r_value,sizeof(int64_t));
    CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);

    // mpc->GMW_A->open(r_value,r_value);


    if (party != 1) {
			io.recv_data(1, &mask_input, 64*sizeof(bool));
			io.flush(1);
		}
		else {
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool.enqueue([&io, mask_input, party2]() {
					io.send_data(party2, &mask_input, 64*sizeof(bool));
					io.flush(party2);
				}));
			}
			joinNclean(res);
		}

    if (party == 1)
    {
        uint64_t band2 = io.count();
	    cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
    }
    


}
