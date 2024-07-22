#include <emp-tool/emp-tool.h>
#include "../../MPABY_CORE/CORE_protocol.h"
#include "../../emp-agmpc/MPGC.h"
#include <omp.h>
// #include "../../MPABY_GC/GC_mpc.h"
using namespace std;
using namespace emp;


#define BENCH 1
#define MMTR_TEST
// #define TR_TEST
// #define MM_TEST
const static int nP = 32;
int party, port;

const static int rows = 10;


const static uint64_t BB = 1;
const static uint64_t feature = 1;
const static uint64_t hidden = 1;
const static uint64_t length = BB*feature;
const static uint64_t length1 = feature*hidden;



// const char out3[] = "7e50ebb4bd718d7b";


int main(int argc, char** argv) {

    Eigen::setNbThreads(omp_get_max_threads());
    std::cout << "Number of threads: " << Eigen::nbThreads() << std::endl;

	parse_party_and_port(argv, &party, &port);
	printf("party:%d	port:%d\n",party,port);
    if(party > nP)return 0; 
    struct timeval t1,t2;
    double timeuse = 0;
    struct timeval total1,total2;
    double timeuse_total = 0;

    NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};

	
	ThreadPool pool(nP);	
    PRG prg;
    

    #ifdef TR_TEST

    int64_t Trunc_A = 0x3000001;
    int64_t Trunc_B = 0;
    int64_t Trunc_C = 0;
    // prg.random_data(&Trunc_A,sizeof(int64_t));
    for (int i = 0; i < 10; i++)
    {
        prg.random_data(&Trunc_B,sizeof(int64_t));
        prg.random_data(&Trunc_A,sizeof(int16_t));

        Trunc_C = Trunc_A + Trunc_B;

        int64_t temp = (Trunc_C >> 13) - ((Trunc_A >> 13) + (Trunc_B >> 13));
        if(abs(temp) == 1 || temp == 0) cout << "success" <<endl;
        else    
        {
            cout << "fail" << endl;
            cout << (Trunc_C >> 13) << "    ";
            cout << (Trunc_A >> 13) + (Trunc_B >> 13) << endl;
        }
        /* code */
    }


    uint64_t a_value[length],b_value[feature];


    prg.random_data(a_value,length*sizeof(uint64_t));
    prg.random_data(b_value,length*sizeof(uint64_t));

    for (int i = 0; i < length; i++)
        {
            a_value[i] = a_value[i] & (0xFFFFFFFFFFFFFFFF>>(64-14));
            /* code */
        }
    // for (int i = 0; i < feature; i++)
    //     {
    //         b_value[i] = b_value[i] & (0xFFFFFFFFFFFFFFFF>>(64-7));
    //         /* code */
    //     }
        

    Map<MatrixX<uint64_t>> A(a_value,BB,feature);
    Map<MatrixX<uint64_t>> B(b_value,BB,feature);
    // MatrixX<int64_t> TEMP = A*B;
    cout << "A:"<<std::hex<<(A/256+B/256)<<endl;
    cout << "B:"<<std::hex<<(A+B)/256<<endl;
    // cout << "C:"<<std::hex<<TEMP/256<<endl;
    #endif
    
	

    
    // NetIOMP<nP> io(party, port);

    //test for vec protocol
    #ifdef VEC
    for (int i = 0; i < 1; i++)
    {   
        	
        EMPDM_protocol<nP>* mpc = new EMPDM_protocol<nP>(&io,&pool,party);
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


        prg.random_data(a_value,length*sizeof(int64_t));
        prg.random_data(b_value,length*sizeof(int64_t));

        for (int i = 0; i < length; i++)
        {
            mpc->GMW_A->set_ShareA(Delta_a[i],delta_a[i],a_value[i],2);
            // delta_a[i] = A[i].delta;
            mpc->GMW_A->set_ShareA(Delta_b[i],delta_b[i],b_value[i],2);
            // delta_b[i] = B[i].delta;
        }

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

        // RowVectorX<int64_t> Vec(Delta_A.cols());
        // Vec.setRandom();
        // cout<<"row:"<<Vec<<endl;
        // cout<<"row:"<<Vec.cols()<<endl;
        // cout<<"row:"<<Vec.col(0)<<endl;

        // Delta_D(0,0) = (RowVectorX<int64_t>::Ones(Delta_A.cols()) * Vec.transpose());
        // cout<<"D:"<<Delta_D<<endl;
        mpc->mpc->MM_setup(Delta_C,Delta_AB,Delta_A,Delta_B);
        mpc->mpc->MM_online(PublicVal_C,Delta_C,Delta_AB,PublicVal_A,PublicVal_B,Delta_A,Delta_B);
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
    }

    #endif

    #ifdef MMTR_TEST

        CORE_protocol<nP>* mpc = new CORE_protocol<nP>(&io,ios,&pool,party);
        // cout <<"Setup MPC:\t"<<party<<"\n";
        int64_t a_value[length],b_value[length1];

        int64_t delta_a[length] = {0};
        int64_t delta_b[length1] = {0};

        int64_t Delta_a[length] = {0};
        int64_t Delta_b[length1] = {0};


        prg.random_data(a_value,length*sizeof(int64_t));
        prg.random_data(b_value,length1*sizeof(int64_t));

        for (uint64_t i = 0; i < length; i++)
        {
            a_value[i] = a_value[i] >> 44;
            /* code */
        }
        for (uint64_t i = 0; i < length1; i++)
        {
            b_value[i] = b_value[i] >> 44;
            /* code */
        }
        

        Map<MatrixX<int64_t>> A(a_value,BB,feature);
        Map<MatrixX<int64_t>> B(b_value,feature,hidden);
        
        // cout << "A:"<<std::hex<<A<<endl;
        // cout << "B:"<<std::hex<<B<<endl;
        // cout << "A*B:"<<std::hex<<A*B<<endl;
        // cout << "C:"<<std::hex<<TEMP<<endl;

        for (uint64_t i = 0; i < length; i++)
        {
            mpc->GMW_A->set_ShareA(Delta_a[i],delta_a[i],a_value[i],1);
            // delta_a[i] = A[i].delta;
            
        }
        for (uint64_t i = 0; i < length1; i++)
        {
            mpc->GMW_A->set_ShareA(Delta_b[i],delta_b[i],b_value[i],2);
            // delta_b[i] = B[i].delta;
        }

        Map<MatrixX<int64_t>> PublicVal_A(Delta_a,BB,feature);
        Map<MatrixX<int64_t>> PublicVal_B(Delta_b,feature,hidden);
        // MatrixX<int64_t> PublicVal_C(2,length/2);
        
        Map<MatrixX<int64_t>> Delta_A(delta_a,BB,feature);
        Map<MatrixX<int64_t>> Delta_B(delta_b,feature,hidden);

        MatrixX<int64_t> Delta_AB(BB,hidden);
        MatrixX<int64_t> Delta_R_tr(BB,hidden);
        MatrixX<int64_t> Delta_R(BB,hidden);
        MatrixX<int64_t> Public_R_tr(BB,hidden);

        MatrixX<int64_t> Public_C(BB,hidden);
        MatrixX<int64_t> Delta_C(BB,hidden);

        MatrixX<int64_t> Public_Rmsb(BB,hidden);
        MatrixX<int64_t> Delta_Rmsb(BB,hidden);

        
        mpc->MPABY_MMTR_setup(Delta_C,Delta_Rmsb,Public_R_tr,Delta_R_tr,Delta_R,Delta_AB,Delta_A,Delta_B,10);
        mpc->MPABY_MMTR_online(Delta_Rmsb,Public_C,Delta_C,Public_R_tr,Delta_R_tr,Delta_R,Delta_AB,PublicVal_A,Delta_A,PublicVal_B,Delta_B,10);

        int64_t c_value[Public_C.size()];
        mpc->GMW_A->rec_ShareA_vec(c_value,Public_C.data(),Delta_C.data(),Public_C.size());

        mpc->GMW_A->rec_ShareA_vec(A.data(),PublicVal_A.data(),Delta_A.data(),PublicVal_A.size());
        mpc->GMW_A->rec_ShareA_vec(B.data(),PublicVal_B.data(),Delta_B.data(),PublicVal_B.size());
        MatrixX<int64_t> TEMP = (A*B).array();
        // cout <<"TEMP"<< endl<<std::hex << TEMP << endl;
        TEMP = (A*B)/((uint64_t)1<<10);
        
        for (int i = 0; i < Public_C.size(); i++)
        {
            int64_t temp = TEMP.data()[i] - c_value[i];
            if(abs(temp) == 1 || temp == 0) cout << "success" <<endl;
            else    
            {
                cout << "fail" << endl;
                cout << std::hex <<TEMP.data()[i] << "    ";
                cout << std::hex <<c_value[i] << endl;
            }
        }
        cout<< endl;


    for (int i = 0; i < BENCH; i++)
    {
        
        
        gettimeofday(&total1,NULL);
        mpc->MPABY_MMTR_setup(Delta_C,Delta_Rmsb,Public_R_tr,Delta_R_tr,Delta_R,Delta_AB,Delta_A,Delta_B,10);
        gettimeofday(&t1,NULL);
        mpc->MPABY_MMTR_online(Delta_Rmsb,Public_C,Delta_C,Public_R_tr,Delta_R_tr,Delta_R,Delta_AB,PublicVal_A,Delta_A,PublicVal_B,Delta_B,10);
        // mpc->GMW_A->rec_ShareA_vec(c_value,Public_C.data(),Delta_C.data(),Public_C.size());
        // mpc->EMPDM_MMTR_online(Public_C,Delta_C,Public_R_tr,Delta_R_tr,Delta_R,Delta_AB,PublicVal_A,Delta_A,PublicVal_B,Delta_B,2);
        gettimeofday(&t2,NULL);
        gettimeofday(&total2,NULL);
        timeuse =  timeuse + (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        timeuse_total =  timeuse_total + (total2.tv_sec - total1.tv_sec) + (double)(total2.tv_usec - total1.tv_usec)/1000000.0;

    }
    // cout << "throughput:"<<(BENCH/timeuse)<<"opt/s"<<endl;
    cout<<"MMTR online runtime: "<<(timeuse/BENCH)*1000.0<<"ms"<<endl;
    cout<<"MMTR_total runtime: "<<(timeuse_total/BENCH)*1000.0<<"ms"<<endl;
    #endif
}


