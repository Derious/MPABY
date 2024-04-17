#include <emp-tool/emp-tool.h>
#include "../../MPABY_GMW/GMW_protocol.h"
using namespace std;
using namespace emp;


const static int nP = 3;
int party, port;

int main(int argc, char** argv) {
 
    int port, party;
	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
	if(party > nP)return 0;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io1(party, port);
	ThreadPool pool(4);	
    PRG prg;
	GMWprotocolA<nP> *gmw = new GMWprotocolA<nP>(&io,&pool,party);
	GMWprotocolB<nP> *gmwb = new GMWprotocolB<nP>(&io1,&pool,party);

	for (int i = 0; i < 10; i++)
	{
		int64_t a_value, b_value, delta_c = 0, delta_ab = 0;
		int64_t rec_value = 0;
		ShareA a, b, c;
		// a_value = 5493878114719452298;
		b_value = 5493878114719452297;
		prg.random_data(&a_value,sizeof(int64_t));
		// prg.random_data(&b_value,sizeof(int64_t));

		gmw->set_ShareA(&a,a_value,1);
		gmw->set_ShareA(&b,b_value,2);

		gmw->multA_setup(delta_c,delta_ab,&a,&b);
		gmw->multA_online(&c,&a,&b,delta_c,delta_ab);
		gmw->rec_ShareA(rec_value,&c);

		cout<<"rec_value:"<<rec_value<<"\tc:"<<a_value*b_value<<endl;

		ShareA r;
		gmw->randBit(&r);
		gmw->rec_ShareA(rec_value,&r);
		cout<<"random_value:"<<rec_value<<endl;
		// gmw->open(rec_value,a.delta);
		// gmw->open(rec_value,b.delta);
		// gmw->addA(&c,&a,&b);
		// gmw->rec_ShareA(rec_value,&c);
		// gmw->setupMULT(rec_value,&a,&b);
		// cout<<"Public Delta:"<<c.PublicVal<<endl;
		// cout<<"Share Delta:"<<c.delta<<endl;
		// cout<<"rec:"<<rec_value<<"\tc:"<<a_value+b_value<<endl;
	}
	
	

	// for (int i = 0; i < 10; i++)
	// {
	// 	bool b_value;
	// 	bool recb_value;
	// 	ShareB b;
	// 	prg.random_bool(&b_value,sizeof(bool));
	// 	gmwb->set_ShareB(&b,b_value,2);
	// 	gmwb->open(recb_value,b.delta);
	// 	// gmwb->rec_ShareB(recb_value,&b);
	// 	// cout<<"recb:"<<recb_value<<"	b:"<<b_value<<endl;
	// }
	
	

}