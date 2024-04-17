#ifndef GMWPROT_H__
#define GMWPROT_H__


#include "../MPABY_util/util_netmp.h"
#include "../MPABY_util/util_helper.h"
#include "../emp-agmpc/MPGC.h"
#include <Eigen/Dense>

using namespace emp;

// #define _debug

typedef struct SHAREA
{
    int64_t PublicVal = 0;
    int64_t delta = 0;  
} ShareA;

typedef struct SHAREB
{
    bool PublicVal;
    bool delta; 
} ShareB;


typedef struct ATRIPLE
{
    int64_t delta_a = 0;
    int64_t delta_b = 0;
    int64_t delta_ab = 0;
} ATriple;

typedef struct A2BCTX
{
    int64_t d = 0;
    bool d_delta[64] = {0};   //input length
    ShareB ri_B[64];    //input length
    ShareA r_A;
    bool res_delta[65] = {0};

} A2B_ctx ;

typedef struct A2BCTX2
{
    // int64_t d = 0;
    // bool d_delta[64] = {0};   //input length
    // ShareB ri_B[64];    //input length
    // ShareA r_A;
    bool input_value[128] = {0};
    bool res_delta[65] = {0};

} A2B_ctx2 ;

typedef struct BIT2ACTX
{
    bool c; 
    ShareA r;

} Bit2A_ctx ;

typedef struct BIT2ACTX2
{
    bool c; 
    int64_t res_delta;
    ShareA r;

} Bit2A_ctx2 ;

typedef struct PQCTX
{
    int64_t delta_c = 0;
    int64_t delta_pq = 0;
    bool c1 = 0;
    ShareA r_p;
    bool c2 = 0;
    ShareA r_q;
    /* data */
} PQ_ctx;

typedef struct PVCTX
{
    int64_t delta_c = 0;
    int64_t delta_pv = 0;
    bool c;
    ShareA r_p;
} PV_ctx;

typedef struct PQVCTX
{
    PV_ctx ctx[2];

} PQV_ctx;



template<int length>
struct scalar_ctx
{
    int64_t delta_c = 0;
    int64_t delta_ab[length] = {0};
};




template<int nP>
class GMWprotocolA
{

public:
    NetIOMP<nP> *io;
	ThreadPool * pool;
	int party;
    PRG prg;
    GMWprotocolA(NetIOMP<nP>* io, ThreadPool * pool, int party){
        this->io = io;
        this->pool = pool;
        this->party = party;
    }

    void set_ShareA(ShareA* a, int64_t value, int p) {

            vector<future<void>> res;//relic multi-thread problems...

            if(p == party) {
                prg.random_data(&(a->delta),8);
                a->PublicVal = a->delta + value;
                for (int i = 1; i <= nP; i++) {
                    int party2 = i;
                    if (i != party) {
                        res.push_back(pool->enqueue([this,a,party2]() {
                            io->send_data(party2,&(a->PublicVal),sizeof(int64_t));
                            io->flush(party2);
				            }));
                    }
                }
                joinNclean(res);
            }
            else {
                a->delta = 0;
                io->recv_data(p,&(a->PublicVal),sizeof(int64_t));
                io->flush(p);
            }

#ifdef _debug_

            cout<<"Public Delta:"<<a->PublicVal<<endl;
            cout<<"Share delta:"<<a->delta<<endl;

#endif

    }

    void set_ShareA(int64_t &PublicVal, int64_t &delta, int64_t value, int p) {

            vector<future<void>> res;//relic multi-thread problems...

            if(p == party) {
                prg.random_data(&(delta),8);
                PublicVal = delta + value;
                for (int i = 1; i <= nP; i++) {
                    int party2 = i;
                    if (i != party) {
                        res.push_back(pool->enqueue([this,PublicVal,party2]() {
                            io->send_data(party2,&(PublicVal),sizeof(int64_t));
                            io->flush(party2);
				            }));
                    }
                }
                joinNclean(res);
            }
            else {
                delta = 0;
                io->recv_data(p,&(PublicVal),sizeof(int64_t));
                io->flush(p);
            }

#ifdef _debug_

            cout<<"Public Delta:"<<a->PublicVal<<endl;
            cout<<"Share delta:"<<a->delta<<endl;

#endif

    }

    void set_ShareA_vec(int64_t* PublicVal, int64_t* delta, int64_t* value, int length, int p) {

            vector<future<void>> res;//relic multi-thread problems...

            if(p == party) {
                prg.random_data(delta,length*sizeof(int64_t));
                for (int i = 0; i < length; i++)
                {
                    PublicVal[i] = delta[i] + value[i];
                }
                for (int i = 1; i <= nP; i++) {
                    int party2 = i;
                    if (i != party) {
                        res.push_back(pool->enqueue([this,PublicVal,party2,length]() {
                            io->send_data(party2,PublicVal,length*sizeof(int64_t));
                            io->flush(party2);
				            }));
                    }
                }
                joinNclean(res);
            }
            else {
                delta = 0;
                io->recv_data(p,PublicVal,length*sizeof(int64_t));
                io->flush(p);
            }

#ifdef _debug_

            cout<<"Public Delta:"<<a->PublicVal<<endl;
            cout<<"Share delta:"<<a->delta<<endl;

#endif

    }

    void rec_ShareA(int64_t &value, ShareA* a){
        
        int64_t tmp = 0;
        int64_t *delta_i = new int64_t[nP+1];
        delta_i[party] = a->delta;

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,a,party2]() {
                    io->send_data(party2,&(a->delta),sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,delta_i,party2]() {
                    io->recv_data(party2,delta_i+party2,sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);
        // cout<<"==================================="<<endl;
        // cout<<"PublicVal:"<<a->PublicVal<<endl;
        for (int i = 1; i <= nP; i++)
        {
            tmp += delta_i[i];
            // cout<<"delta_i:"<<delta_i[i]<<" tmp:"<<tmp<<endl;
        }
        // cout<<"==================================="<<endl;
        value = a->PublicVal - tmp;
        // cout<<"value:"<<value<<endl;

        delete[] delta_i;
        
#ifdef _debug
        cout<<"value:"<<value<<endl;
#endif
    }

    void rec_ShareA_vec(int64_t* value, int64_t* Publicvalue, int64_t* share_value, int length) {

        int64_t *tmp[nP+1] = {nullptr};

        for (int i = 1; i <= nP; i++)
        {
            tmp[i] = new int64_t[length];
            /* code */
        }
        
        int64_t sum[length] = {0};

        memcpy(tmp[party],share_value,length*sizeof(int64_t));

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,share_value,length,party2]() {
                    io->send_data(party2,share_value,length*sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,tmp,length,party2]() {
                    io->recv_data(party2,tmp[party2],length*sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int j = 0; j < length; j++)
        {
            for (int i = 1; i <= nP; i++)
            {
                sum[j] += tmp[i][j];
                // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
                #ifdef _debug_
                cout<<"share_value:"<<tmp[i]<<endl;
                #endif
            }
            #ifdef _debug_
            cout<<"sum:"<<sum<<endl;
            #endif
            value[j] = Publicvalue[j] - sum[j];
        }
    }


    //todo: add sp19 randBit protocol
    void randBit(ShareA* r) {

        bool random_value = 0;

        if (party == 1) {
            prg.random_bool(&random_value,1);
        }

        set_ShareA(r,(int64_t)random_value,1);
        
    }

    //todo: add sp19 randBit protocol
    void randBit(int64_t* PublicVal, int64_t* delta, int length) {

        bool random_value[length] = {0};

        if (party == 1) {
            prg.random_bool(random_value,length);
        }
        for (int i = 0; i < length; i++)
        {   
            // cout<<i<<":"<<random_value[i]<<endl;
            set_ShareA(PublicVal[i], delta[i], (int64_t)random_value[i],1);
            /* code */
        }   
    }

    void open(int64_t &value, int64_t share_value) {

        int64_t* tmp = new int64_t[nP+1];
        int64_t sum = 0;
        tmp[party] = share_value;

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,share_value,party2]() {
                    io->send_data(party2,&share_value,sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,tmp,party2]() {
                    io->recv_data(party2,tmp+party2,sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int i = 1; i <= nP; i++)
        {
            sum += tmp[i];
            // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
            #ifdef _debug_
            cout<<"share_value:"<<tmp[i]<<endl;
            #endif
        }
        #ifdef _debug_
            cout<<"sum:"<<sum<<endl;
        #endif
        value = sum;

    }

    void open_vec(int64_t* value, int64_t* share_value, int length) {

        int64_t *tmp[nP+1] = {nullptr};

        for (int i = 1; i <= nP; i++)
        {
            tmp[i] = new int64_t[length];
            /* code */
        }
        
        int64_t sum[length] = {0};

        memcpy(tmp[party],share_value,length*sizeof(int64_t));

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,share_value,length,party2]() {
                    io->send_data(party2,share_value,length*sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,tmp,length,party2]() {
                    io->recv_data(party2,tmp[party2],length*sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int j = 0; j < length; j++)
        {
            for (int i = 1; i <= nP; i++)
            {
                sum[j] += tmp[i][j];
                // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
                #ifdef _debug_
                cout<<"share_value:"<<tmp[i]<<endl;
                #endif
            }
            #ifdef _debug_
            cout<<"sum:"<<sum<<endl;
            #endif
            value[j] = sum[j];
        }
    }

    void addA_setup(int64_t &delta_c, int64_t delta_a, int64_t delta_b) {
        delta_c = delta_a + delta_b;
    }

    void addA_online(ShareA*c, ShareA* a, ShareA* b) {
        c->PublicVal = a->PublicVal + b->PublicVal;
    }

    void addA_setup_vec(Eigen::RowVectorX<int64_t> &delta_c, Eigen::RowVectorX<int64_t> delta_a, Eigen::RowVectorX<int64_t> delta_b) {
        delta_c = delta_a + delta_b;
    }
    void addA_online_vec(Eigen::RowVectorX<int64_t> &Delta_c, Eigen::RowVectorX<int64_t> Delta_a, Eigen::RowVectorX<int64_t> Delta_b) {
        Delta_c = Delta_a + Delta_b;
    }

    void addA_setup_mat(Eigen::MatrixX<int64_t> &delta_c, Eigen::MatrixX<int64_t> delta_a, Eigen::MatrixX<int64_t> delta_b) {
        delta_c = delta_a + delta_b;
    }
    void addA_online_mat(Eigen::MatrixX<int64_t> &Delta_c, Eigen::MatrixX<int64_t> Delta_a, Eigen::MatrixX<int64_t> Delta_b) {
        Delta_c = Delta_a + Delta_b;
    }

    void addA(ShareA*c, ShareA* a, ShareA* b) {

        c->delta = a->delta + b->delta;
        c->PublicVal = a->PublicVal + b->PublicVal;

    }

    void addA_constant_setup(int64_t &delta_c, int64_t delta_a, int64_t b) {
        delta_c = delta_a ;
    }

    void addA_constant_online(ShareA*c, ShareA* a, int64_t b) {
        c->PublicVal = a->PublicVal + b;
    }

    void addA_constant_setup_vec(Eigen::RowVectorX<int64_t> &delta_c, Eigen::RowVectorX<int64_t> delta_a, Eigen::RowVectorX<int64_t> b) {
        delta_c = delta_a ;
    }
    void addA_constant_online_vec(Eigen::RowVectorX<int64_t> &Delta_c, Eigen::RowVectorX<int64_t> Delta_a, Eigen::RowVectorX<int64_t> b) {
        Delta_c = Delta_a + b;
    }

    void addA_constant(ShareA*c, ShareA* a, int64_t b) {

            c->delta = a->delta;
            c->PublicVal = a->PublicVal + b;
    }

    void multA_setup(int64_t &delta_c, int64_t &delta_ab, ShareA* a, ShareA* b) {

        prg.random_data(&delta_c,sizeof(int64_t)); // generate delta_c

        setupMULT(delta_ab,a,b);
    }

    void multA_setup(int64_t &delta_c, int64_t &delta_ab, int64_t &a_delta, int64_t &b_delta) {

        prg.random_data(&delta_c,sizeof(int64_t)); // generate delta_c

        setupMULT(delta_ab,a_delta,b_delta);
    }

    void multA_online(ShareA* c, ShareA* a, ShareA* b, int64_t delta_c, int64_t delta_ab) {

        int64_t Delta_share = 0;
        int64_t Delta_sum = 0;
        if (party == 1)
            Delta_share = a->PublicVal*b->PublicVal - a->PublicVal*b->delta - b->PublicVal*a->delta + delta_ab +delta_c;
        else
            Delta_share = delta_ab +delta_c - a->PublicVal*b->delta - b->PublicVal*a->delta;

        
        open(Delta_sum,Delta_share);

        c->delta = delta_c;
        c->PublicVal = Delta_sum;     
    }


    /*
    Without considering the implementation of OT, calculate it as a black box. Currently, 
    there is no need to consider the communication volume of setup. 
    This is just a simple implementation.
    */
    void setupMULT(int64_t &delta_ab, ShareA* a, ShareA* b) {

        
        int64_t *delta_i = new int64_t[nP+1];
        int64_t *delta2_i = new int64_t[nP+1];

        int64_t deltaa = 0,deltab = 0,temp = 0;
        delta_i[party] = a->delta;
        delta2_i[party] = b->delta;


        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,a,b,party2]() {
                    io->send_data(party2,&(a->delta),sizeof(int64_t));
                    io->flush(party2);
                    io->send_data(party2,&(b->delta),sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,delta_i,delta2_i,party2]() {
                    io->recv_data(party2,delta_i+party2,sizeof(int64_t));
                    io->flush(party2);
                    io->recv_data(party2,delta2_i+party2,sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int i = 1; i <= nP; i++)
        {
            deltaa += delta_i[i];
            deltab += delta2_i[i];
            // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
        }
        temp = deltaa*deltab;

        if (party == 1)
            delta_ab = temp;
        else
            delta_ab = 0;
    }

    void setupMULT(int64_t &delta_ab, int64_t &a_delta, int64_t &b_delta) {

        
        int64_t *delta_i = new int64_t[nP+1];
        int64_t *delta2_i = new int64_t[nP+1];

        int64_t deltaa = 0,deltab = 0,temp = 0;
        delta_i[party] = a_delta;
        delta2_i[party] = b_delta;


        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,a_delta,b_delta,party2]() {
                    io->send_data(party2,&(a_delta),sizeof(int64_t));
                    io->flush(party2);
                    io->send_data(party2,&(b_delta),sizeof(int64_t));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,delta_i,delta2_i,party2]() {
                    io->recv_data(party2,delta_i+party2,sizeof(int64_t));
                    io->flush(party2);
                    io->recv_data(party2,delta2_i+party2,sizeof(int64_t));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int i = 1; i <= nP; i++)
        {
            deltaa += delta_i[i];
            deltab += delta2_i[i];
            // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
        }
        temp = deltaa*deltab;

        if (party == 1)
            delta_ab = temp;
        else
            delta_ab = 0;
    }

    void mulA_constant_setup(int64_t &delta_c, int64_t delta_a, int64_t b) {
        delta_c = b * delta_a;
    }

    void mulA_constant_online(ShareA*c, ShareA* a, int64_t b) {
        c->PublicVal = b * a->PublicVal;
    }


    void mulA_constant_setup_vec(Eigen::RowVectorX<int64_t> &delta_c, Eigen::RowVectorX<int64_t> delta_a, int64_t b) {
        delta_c = b * delta_a ;
    }
    void mulA_constant_online_vec(Eigen::RowVectorX<int64_t> &Delta_c, Eigen::RowVectorX<int64_t> Delta_a, int64_t b) {
        Delta_c = b * Delta_a ;
    }

    void mulA_constant_setup_mat(Eigen::MatrixX<int64_t> &delta_c, Eigen::MatrixX<int64_t> delta_a, int64_t b) {
        delta_c = b * delta_a ;
    }
    void mulA_constant_online_mat(Eigen::MatrixX<int64_t> &Delta_c, Eigen::MatrixX<int64_t> Delta_a, int64_t b) {
        Delta_c = b * Delta_a ;
    }

    void mulA_constant(ShareA*c, ShareA* a, int64_t b) {
        c->delta = b * a->delta;
        c->PublicVal = b * a->PublicVal;
    }

    ~GMWprotocolA(){

    }
};

template<int nP>
class GMWprotocolB
{
public:
    NetIOMP<nP> *io;
	ThreadPool * pool;
	int party;
    PRG prg;
    GMWprotocolB(NetIOMP<nP>* io, ThreadPool * pool, int party) {
        this->io = io;
        this->pool = pool;
        this->party = party;
    }

    void set_ShareB(ShareB* a, bool value, int p,  bool * _delta = nullptr) {

            vector<future<void>> res;//relic multi-thread problems...

            if(p == party) {
                if (_delta == nullptr)
                    prg.random_bool(&(a->delta),1);
                else
                    a->delta = *_delta;
                a->PublicVal = a->delta != value;
                for (int i = 1; i <= nP; i++) {
                    int party2 = i;
                    if (i != party) {
                        res.push_back(pool->enqueue([this, a, party2]() {
                            io->send_data(party2,&(a->PublicVal),sizeof(bool));
                            io->flush(party2);
				            }));
                    }
                }
                joinNclean(res);
            }
            else {
                a->delta = 0;
                io->recv_data(p,&(a->PublicVal),sizeof(bool));
                io->flush(p);
            }

#ifdef _debug

            cout<<"Public Delta:"<<a->PublicVal<<endl;
            cout<<"Share delta:"<<a->delta<<endl;

#endif     
    }
    void rec_ShareB(bool &value, ShareB* a){
        
        bool tmp = 0;
        bool *delta_i = new bool[nP+1];
        delta_i[party] = a->delta;

        vector<future<void>> Queue;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                Queue.push_back(pool->enqueue([this, a, party2]() {
                    io->send_data(party2,&a->delta,sizeof(bool));
                    io->flush(party2);
				}));
                Queue.push_back(pool->enqueue([this,delta_i,party2]() {
                    io->recv_data(party2,delta_i+party2,sizeof(bool));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(Queue);
        // delete &Queue;

        for (int i = 1; i <= nP; i++)
        {
            tmp = tmp != delta_i[i];
        }
        value = a->PublicVal != tmp;

        delete[] delta_i;
        
#ifdef _debug
        cout<<"value:"<<value<<endl;
#endif
    }

    void rec_ShareB(bool *value, ShareB* a, int length){
        
        
        bool *delta_i[nP+1] = {nullptr};

        for (int i = 1; i <= nP; i++)
        {
            delta_i[i] = new bool[length];
            /* code */
        }
        

        bool *delta = new bool[length];

        for (int i = 0; i < length; i++)
        {
            delta[i] = a[i].delta;
        }
        memcpy(delta_i[party],delta,length*sizeof(bool));
        

        vector<future<void>> qq;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                qq.push_back(pool->enqueue([this, delta, length, party2]() {
                    io->send_data(party2,delta,length*sizeof(bool));
                    io->flush(party2);
				}));
                qq.push_back(pool->enqueue([this,delta_i, length, party2]() {
                    io->recv_data(party2,delta_i[party2],length*sizeof(bool));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(qq);
        // delete &Queue;
        for (int j = 0; j < length; j++)
        {
            bool tmp = 0;
            for (int i = 1; i <= nP; i++)
            {
                tmp = tmp != delta_i[i][j];
            }
            value[j] =  a[j].PublicVal != tmp;
            /* code */
        }

        delete[] delta;
        delta = nullptr;
        for (int i = 1; i < nP+1; i++)
        {
             delete[] delta_i[i];
             delta_i[i] = nullptr;
            /* code */
        }
#ifdef _debug
        cout<<"value:"<<value<<endl;
#endif
    }

    void open(bool &value, bool share_value) {

        bool* tmp = new bool[nP+1];
        bool sum = 0;
        tmp[party] = share_value;

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,share_value,party2]() {
                    io->send_data(party2,&share_value,sizeof(bool));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,tmp,party2]() {
                    io->recv_data(party2,tmp+party2,sizeof(bool));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int i = 1; i <= nP; i++)
        {
            sum = sum != tmp[i];
            // cout<<"delta_a:"<<delta_i[i]<<"\tdelta_b:"<<delta2_i[i]<<"\tdeltab:"<<deltab<<endl;
            #ifdef _debug
            cout<<"share_value:"<<tmp[i]<<endl;
            #endif
        }
        #ifdef _debug
            cout<<"sum:"<<sum<<endl;
        #endif
        value = sum;

    }

    void open_vec(bool* value, bool* share_value, int length) {

        bool *tmp[nP+1] = {nullptr};

        for (int i = 1; i <= nP; i++)
        {
            tmp[i] = new bool[length];
            /* code */
        }
        
        bool sum[length] = {0};

        memcpy(tmp[party],share_value,length*sizeof(bool));

        vector<future<void>> res;//relic multi-thread problems...
        for (int i = 1; i <= nP; i++)
        {   
            int party2 = i;
            if (i!=party)
            {
                res.push_back(pool->enqueue([this,share_value,length,party2]() {
                    io->send_data(party2,share_value,length*sizeof(bool));
                    io->flush(party2);
				}));
                res.push_back(pool->enqueue([this,tmp,length,party2]() {
                    io->recv_data(party2,tmp[party2],length*sizeof(bool));
                    io->flush(party2);
                }));
            }
        }
        joinNclean(res);

        for (int j = 0; j < length; j++)
        {
            for (int i = 1; i <= nP; i++)
            {
                sum[j] = sum[j] != tmp[i][j];

            }
            value[j] = sum[j];
        }
    }

    ~GMWprotocolB() {

    }
};


#endif