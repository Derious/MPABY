#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"
// #include "../../MPABY_GC/GC_mpc.h"
using namespace std;
using namespace emp;


const static int nP = 3;
int party, port;
const static int length = 256;
const static int rows = 16;

// const char out3[] = "7e50ebb4bd718d7b";

void test(Eigen::MatrixX<int64_t> &PublicVal_C) {

    MatrixX<int64_t> Dd_value(2,length/2);
    Dd_value.setRandom();
    PublicVal_C = Dd_value;
}

int main(int argc, char** argv) {
    
    int64_t temp = 10000;

    printf("%.5f\n",(float)temp/(1<<5));

	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
	if(party > nP)return 0;

	NetIOMP<nP> io(party, port);
	ThreadPool pool(4);	
    PRG prg;

    NetIOMP<nP> io_1(party, port+2*(nP+1)*(nP+1)+2);
	NetIOMP<nP> io_2(party, port+2*(nP+1)*(nP+1)+3);
	NetIOMP<nP> *ios[2] = {&io_1, &io_2};


    //test for vec protocol
    for (int i = 0; i < 1; i++)
    {   
        	
        CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);
        // cout <<"Setup MPC:\t"<<party<<"\n";

        ShareA res;

        int64_t a_value[length],b_value[length];
        ShareA A[length], B[length];

        int64_t delta_c = 0;
        int64_t delta_ab[length] = {0};
        int64_t delta_a[length] = {0};
        int64_t delta_b[length] = {0};

        int64_t Delta_a[length] = {0};
        int64_t Delta_b[length] = {0};


        // // int64_t delta_c = 0, delta_pq = 0;
        // // bool c1;
        // prg.random_bool(&p_v,1);
        // prg.random_bool(&q_v,1);
        prg.random_data(a_value,length*sizeof(int64_t));
        prg.random_data(b_value,length*sizeof(int64_t));

        for (int i = 0; i < length; i++)
        {
            mpc->GMW_A->set_ShareA(Delta_a[i],delta_a[i],a_value[i],2);
            // delta_a[i] = A[i].delta;
            mpc->GMW_A->set_ShareA(Delta_b[i],delta_b[i],b_value[i],2);
            // delta_b[i] = B[i].delta;
        }
        // Map<RowVectorX<int64_t>> PublicVal_A(Delta_a,length);
        // Map<RowVectorX<int64_t>> PublicVal_B(Delta_b,length);
        // RowVectorX<int64_t> PublicVal_C(length);
        
        // Map<RowVectorX<int64_t>> Delta_A(delta_a,length);
        // Map<RowVectorX<int64_t>> Delta_B(delta_b,length);
        // RowVectorX<int64_t> Delta_C(length);

        Map<MatrixX<int64_t>> PublicVal_A(Delta_a,rows,length/rows);
        Map<MatrixX<int64_t>> PublicVal_B(Delta_b,rows,length/rows);
        // MatrixX<int64_t> PublicVal_C(2,length/2);
        
        Map<MatrixX<int64_t>> Delta_A(delta_a,rows,length/rows);
        Map<MatrixX<int64_t>> Delta_B(delta_b,rows,length/rows);
        // MatrixX<int64_t> Delta_C(2,length/2);
        // MatrixX<int64_t> Delta_D(2,length/2);

        MatrixX<int64_t> Delta_C(rows,length/rows);
        MatrixX<int64_t> Delta_AB(rows,length/rows);
        MatrixX<int64_t> PublicVal_C(rows,length/rows);

        srand(clock());

        // for (int i = 0; i < 100; i++)
        // {
        //     Delta_D.setRandom();
        //     cout<<"D:"<< std::hex << Delta_D<<endl;
        //     /* code */
        // }


        // for (int i = 0; i < Delta_A.rows(); i++)
        // {
        //     for (int j = 0; j < Delta_B.cols(); j++)
        //     {
        //         for (int k = 0; k < Delta_A.cols(); k++)
        //         {
        //             cout<<"("<<i<<" , "<<k<<")";
        //             cout<<" ("<<k<<" , "<<j<<")"<<endl;
        //         }
        //         cout<<"location: ("<<i<<" , "<<j<<")"<<endl;
        //     }
        // }
        // // cout<<"row:"<<Delta_D.rows()<<endl;
        // // cout<<"row:"<<Delta_D.cols()<<endl;
        // // cout<<"D:"<<Delta_D(0,0)<< std::hex <<endl;

        // RowVectorX<int64_t> Vec(Delta_A.cols());
        // Vec.setRandom();
        // cout<<"row:"<<Vec<<endl;
        // cout<<"row:"<<Vec.cols()<<endl;
        // cout<<"row:"<<Vec.col(0)<<endl;

        // Delta_D(0,0) = (RowVectorX<int64_t>::Ones(Delta_A.cols()) * Vec.transpose());
        // cout<<"D:"<<Delta_D<<endl;
        mpc->MM_setup(Delta_C,Delta_AB,Delta_A,Delta_B);
        mpc->MM_online(PublicVal_C,Delta_C,Delta_AB,PublicVal_A,PublicVal_B,Delta_A,Delta_B);
        // cout<<"PublicVal_C:"<<PublicVal_C<<endl;
        // cout<<"Delta_C:"<<Delta_C<<endl;
        int64_t Cc_rec[length] = {0};
        mpc->GMW_A->rec_ShareA_vec(Cc_rec, PublicVal_C.data(), Delta_C.data(), PublicVal_C.size());
        


        Map<MatrixX<int64_t>> Aa_value(a_value,rows,length/rows);
        Map<MatrixX<int64_t>> Bb_value(b_value,rows,length/rows);
        MatrixX<int64_t> Cc_value(rows,length/rows);
        Cc_value = Aa_value*Bb_value;
        cout<<"Aa_value:"<<Aa_value<<endl;
        cout<<"Bb_value:"<<Bb_value<<endl;
        cout<<"Cc_value:"<<Cc_value<<endl;

        for (int i = 0; i < length; i++)
        {
            cout<<i<<":"<<(Cc_rec[i] == Cc_value.data()[i]?"success":"fail")<<endl;
        }


        // int64_t aBit[length] = {0};
        // int64_t Public_aBit[length] = {0};
        // int64_t Delta_aBit[length] = {0};
        // int64_t aBit_rec[length] = {0};
        // mpc->GMW_A->randBit(Public_aBit,Delta_aBit,length);

        // mpc->GMW_A->rec_ShareA_vec(aBit_rec,Public_aBit, Delta_aBit, length);
        // for (int i = 0; i < length; i++)
        // {
        //     cout<<i<<":"<<aBit_rec[i]<<endl;
        // }
        


        // MatrixX<int64_t> Dd_value(2,length/2);
        // test(Dd_value);
        // cout<<"Dd_value:"<<Dd_value<<endl;



    }
    

}

