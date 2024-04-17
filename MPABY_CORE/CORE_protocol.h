#ifndef MPABYPROT_H__
#define MPABYPROT_H__

#define _debug_core
#include "../MPABY_GMW/GMW_protocol.h"
#include "../emp-agmpc/MPGC.h"
#include <omp.h>
// #include <Eigen/Dense>
#include <sys/time.h>
// #include "../MPABY_GC/GC_mpc.h"
#define LSB(a) (uint64_t)(a) & (uint64_t)0x1
using namespace emp;
using namespace Eigen;


template<int nP>
class CORE_protocol
{
public:

    NetIOMP<nP> *io;
	ThreadPool * pool;
	int party;
    PRG prg;
    GMWprotocolA<nP> *GMW_A;
    GMWprotocolB<nP> *GMW_B;
    NetIOMP<nP> *ios[2];
    CORE_protocol(NetIOMP<nP>* io, NetIOMP<nP> * IOs[2], ThreadPool * pool, int party){
        this->io = io;
        this->pool = pool;
        this->party = party;
        // this->ios[0] = IOs[0];
        // this->ios[1] = IOs[1];
        GMW_A = new GMWprotocolA<nP>(io,pool,party);
        GMW_B = new GMWprotocolB<nP>(io,pool,party);
    }

    //no interaction， communication-free
    void A2Bit(ShareB* res, ShareA* a) {
            res->delta = LSB(a->delta);
            res->PublicVal = LSB(a->PublicVal);
    }

    void Bit2A_setup(bool &c, ShareA* r, ShareB* a) {


        bool Ri, Ci;
        GMW_A->randBit(r);

        // cout << "delta:"<<r->delta<<"\tpublic:"<<r->PublicVal<<endl;


        if (party == 1) 
            Ri = LSB(r->delta)^LSB(r->PublicVal);
        else
            Ri = LSB(r->delta);

        Ci = Ri != a->delta;
        // cout<<"Ci:"<<Ci<<endl;

        GMW_B->open(c,Ci);

        #ifdef _debug_core_
        cout<<"==================================="<<endl;
            cout<<"Bit2A_setup:"<<endl;
            int64_t rec_r = 0;
            bool rec_a = 0;
            GMW_A->rec_ShareA(rec_r,r);
            GMW_B->open(rec_a,a->delta);
            cout << "c:"<<c<<"\tr:"<<rec_r<<"\ta:"<<rec_a<<endl;
            cout <<(c != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
        cout<<"==================================="<<endl;
        #endif
    }

    void Bit2A_setup(bool &c, ShareA* r, bool a_delta) {


        bool Ri, Ci;
        GMW_A->randBit(r);

        // cout << "delta:"<<r->delta<<"\tpublic:"<<r->PublicVal<<endl;


        if (party == 1) 
            Ri = LSB(r->delta)^LSB(r->PublicVal);
        else
            Ri = LSB(r->delta);

        Ci = Ri != a_delta;
        // cout<<"Ci:"<<Ci<<endl;

        GMW_B->open(c,Ci);

        #ifdef _debug_core_
        cout<<"==================================="<<endl;
            cout<<"Bit2A_setup:"<<endl;
            int64_t rec_r = 0;
            bool rec_a = 0;
            GMW_A->rec_ShareA(rec_r,r);
            GMW_B->open(rec_a,a_delta);
            cout << "c:"<<c<<"\tr:"<<rec_r<<"\ta:"<<rec_a<<endl;
            cout <<(c != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
        cout<<"==================================="<<endl;
        #endif
    }
    
    void Bit2A_setup(bool &c, ShareA* r, int64_t &a_deltaA, bool a_delta) {


        bool Ri, Ci;
        GMW_A->randBit(r);
        prg.random_data(&a_deltaA,sizeof(int64_t));

        // cout << "delta:"<<r->delta<<"\tpublic:"<<r->PublicVal<<endl;


        if (party == 1) 
            Ri = LSB(r->delta)^LSB(r->PublicVal);
        else
            Ri = LSB(r->delta);

        Ci = Ri != a_delta;
        // cout<<"Ci:"<<Ci<<endl;

        GMW_B->open(c,Ci);

        #ifdef _debug_core_
        cout<<"==================================="<<endl;
            cout<<"Bit2A_setup:"<<endl;
            int64_t rec_r = 0;
            bool rec_a = 0;
            GMW_A->rec_ShareA(rec_r,r);
            GMW_B->open(rec_a,a_delta);
            cout << "c:"<<c<<"\tr:"<<rec_r<<"\ta:"<<rec_a<<endl;
            cout <<(c != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
        cout<<"==================================="<<endl;
        #endif
    }

    void Bit2A_setup_parallel(bool* c, ShareA* r, int64_t* a_deltaA, bool* a_delta, int num) {


        bool Ri;
        bool Ci[num];
        for (int i = 0; i < num; i++)
        {
            GMW_A->randBit(r);
        }
        
        
        prg.random_data(a_deltaA,num*sizeof(int64_t));

        // cout << "delta:"<<r->delta<<"\tpublic:"<<r->PublicVal<<endl;


        if (party == 1) 
            Ri = LSB(r->delta)^LSB(r->PublicVal);
        else
            Ri = LSB(r->delta);

        for (int i = 0; i < num; i++)
        {
            Ci[i] = Ri != a_delta[i];
        }
        
        
        // cout<<"Ci:"<<Ci<<endl;
        GMW_B->open_vec(c,Ci,num);
        // GMW_B->open(c,Ci);

        #ifdef _debug_core_
        cout<<"==================================="<<endl;
            cout<<"Bit2A_setup:"<<endl;
            int64_t rec_r = 0;
            bool rec_a = 0;
            GMW_A->rec_ShareA(rec_r,r);
            GMW_B->open(rec_a,a_delta);
            cout << "c:"<<c<<"\tr:"<<rec_r<<"\ta:"<<rec_a<<endl;
            cout <<(c != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
        cout<<"==================================="<<endl;
        #endif
    }

    void Bit2A_online(ShareA* res, ShareB* a, bool c, ShareA* r) {
        bool C_prime;
        C_prime = c != a->PublicVal;

        #ifdef _debug_core_
            cout<<"==================================="<<endl;
            cout<<"Bit2A_online:"<<endl;
            int64_t rec_r;
            bool rec_a;
            GMW_B->rec_ShareB(rec_a,a);
            GMW_A->rec_ShareA(rec_r,r);
            // GMW_B->open(rec_a,a->delta);
            cout <<(C_prime != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
            cout << "r:"<<rec_r<<endl;
        #endif

        GMW_A->mulA_constant(res,r,-2);
        #ifdef _debug_core_
            GMW_A->rec_ShareA(rec_r,res);
            cout << "-2r:"<<rec_r<<endl;
        #endif
        GMW_A->addA_constant(res,res,1);
        #ifdef _debug_core_
            GMW_A->rec_ShareA(rec_r,res);
            cout << "-2r+1:"<<rec_r<<endl;
        #endif
        GMW_A->mulA_constant(res,res,(int64_t)C_prime);
        #ifdef _debug_core_
            GMW_A->rec_ShareA(rec_r,res);
            cout << "c(-2r+1):"<<rec_r<<endl;
        #endif
        GMW_A->addA(res,res,r);
        #ifdef _debug_core_
            GMW_A->rec_ShareA(rec_r,res);
            cout << "c(-2r+1)+r:"<<rec_r<<endl;
            cout<<"==================================="<<endl;
        #endif
    }

    void Bit2A_online(ShareA* res, ShareB* a, int64_t res_delta, bool c, ShareA* r) {
        bool C_prime;
        C_prime = c != a->PublicVal;

        #ifdef _debug_core_
            cout<<"==================================="<<endl;
            cout<<"Bit2A_online:"<<endl;
            int64_t rec_r;
            bool rec_a;
            GMW_B->rec_ShareB(rec_a,a);
            GMW_A->rec_ShareA(rec_r,r);
            // GMW_B->open(rec_a,a->delta);
            cout <<(C_prime != (rec_r^(int64_t)rec_a)? "fail":"success")<<endl;
            cout << "r:"<<rec_r<<endl;
        #endif

        //set [[r]]
        int64_t shared_r = 0;
        int64_t shared_res = 0;
        if (party == 1)
            shared_r = r->PublicVal - r->delta;
        else
            shared_r = - r->delta;
        
        shared_res = -2*shared_r;

        if (party == 1)
            shared_res = shared_res + 1;
        else
            shared_r = shared_res;

        shared_res = (int64_t)C_prime*shared_res;

        shared_res = shared_res + shared_r + res_delta;

        res->PublicVal = shared_res;
        // GMW_A->open(res->PublicVal,shared_res);
        res->delta = res_delta;

        // GMW_A->mulA_constant(res,r,-2);
        // #ifdef _debug_core_
        //     GMW_A->rec_ShareA(rec_r,res);
        //     cout << "-2r:"<<rec_r<<endl;
        // #endif
        // GMW_A->addA_constant(res,res,1);
        // #ifdef _debug_core_
        //     GMW_A->rec_ShareA(rec_r,res);
        //     cout << "-2r+1:"<<rec_r<<endl;
        // #endif
        // GMW_A->mulA_constant(res,res,(int64_t)C_prime);
        // #ifdef _debug_core_
        //     GMW_A->rec_ShareA(rec_r,res);
        //     cout << "c(-2r+1):"<<rec_r<<endl;
        // #endif
        // GMW_A->addA(res,res,r);
        // #ifdef _debug_core_
        //     GMW_A->rec_ShareA(rec_r,res);
        //     cout << "c(-2r+1)+r:"<<rec_r<<endl;
        //     cout<<"==================================="<<endl;
        // #endif
    }

    void Bit2A_setup_api(Bit2A_ctx *ctx, bool a_delta) {
        Bit2A_setup(ctx->c, &(ctx->r), a_delta);
    }

    void Bit2A_online_api(ShareA* res, Bit2A_ctx* ctx, ShareB* a) {
        Bit2A_online(res, a, ctx->c, &(ctx->r));
    }

    void Bit2A_setup_api(Bit2A_ctx2 *ctx, bool a_delta) {
        Bit2A_setup(ctx->c, &(ctx->r), ctx->res_delta, a_delta);
    }

    void Bit2A_online_api(ShareA* res, Bit2A_ctx2* ctx, ShareB* a) {
        Bit2A_online(res, a,ctx->res_delta, ctx->c, &(ctx->r));
        GMW_A->open(res->PublicVal,res->PublicVal);
    }

    void B2A_setup_api(Bit2A_ctx ctx[64], bool a_delta[64]) {
        for (int i = 0; i < 64; i++)
        {
            Bit2A_setup_api(ctx+i,a_delta[i]);
        }
    }

    void B2A_setup_api(Bit2A_ctx2 ctx[64],int64_t &res_delta, bool a_delta[64]) {
        int64_t bais = 1;
        int64_t res_deltaA[64] = {0};
        bool C[64] = {0};
        ShareA r;
        Bit2A_setup_parallel(C,&r,res_deltaA,a_delta,64);

        for (int i = 0; i < 64; i++)
        {
            // Bit2A_setup_api(ctx+i,a_delta[i]);
            ctx[i].c = C[i];
            ctx[i].r = r;
            ctx[i].res_delta = res_deltaA[i];
            res_delta = res_delta + bais*res_deltaA[i];
            bais *= 2;
        }
        // for (int i = 0; i < 64; i++)
        // {
        //     Bit2A_setup_api(ctx+i,a_delta[i]);
        //     res_delta = res_delta + bais*((ctx+i)->res_delta);
        //     bais *= 2;
        // }
    }

    void B2A_online_api(ShareA* res, Bit2A_ctx ctx[64], ShareB a[64]) {

        int64_t bais = 1;
        for (int i = 0; i < 64; i++) {
                ShareA tmp;
                Bit2A_online_api(&tmp,ctx+i,a+i);
                GMW_A->mulA_constant(&tmp,&tmp,bais);
                GMW_A->addA(res,res,&tmp);
                bais *= 2;
            }

            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,res);
                printf("rA:%lx\n",rec_A);
            #endif

    }

    void B2A_online_api(ShareA* res, Bit2A_ctx2 ctx[64],int64_t res_delta, ShareB a[64]) {

        int64_t bais = 1;
        int64_t res_public = 0;
        for (int i = 0; i < 64; i++) {
                ShareA tmp;
                // Bit2A_online_api(&tmp,ctx+i,a+i);
                Bit2A_online(&tmp, a+i,(ctx+i)->res_delta, (ctx+i)->c, &((ctx+i)->r));
                res_public = res_public + bais*(tmp.PublicVal);
                // GMW_A->mulA_constant(&tmp,&tmp,bais);
                // GMW_A->addA(res,res,&tmp);
                bais *= 2;
            }
        
        GMW_A->open(res->PublicVal,res_public);
        res->delta = res_delta;

            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,res);
                printf("rA:%lx\n",rec_A);
            #endif

    }


    void A2B_setup(int64_t &d, bool* d_delta, ShareB ri_B[64], ShareA* r_A, ShareA *a, MPGC<nP>* GC) {

            ShareA ri[64];
            int64_t bais = 1;

            // int64_t di[nP+1] = {0};

            int num_in = (GC->cf->n1)+(GC->cf->n2);
            bool* delta = new bool[num_in];
            int64_t* di = new int64_t[nP+1];
            memset(di,0,sizeof(int64_t)*(nP+1));

            for (int i = 0; i < 64; i++) {
                ShareA tmp;
                GMW_A->randBit(ri+i);
                A2Bit(ri_B+i,ri+i);
                GMW_A->mulA_constant(&tmp,ri+i,bais);
                GMW_A->addA(r_A,r_A,&tmp);
                bais *= 2;
            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,ri+i);
                bool rec_B;
                GMW_B->rec_ShareB(rec_B,ri_B+i);
                cout<<(rec_A != (int64_t)rec_B ? "fail" : "success")<<endl;
            #endif
            }


            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,r_A);
                printf("rA:%lx\n",rec_A);
            #endif

            di[party] = a->delta - r_A->delta;

            if (party != 1) {
                io->send_data(1,di+party,sizeof(int64_t));
                io->flush(1);
            }
            else {
                vector<future<void>> res;//relic multi-thread problems...
                for (int i = 2; i <= nP; i++)
                {
                    int party2 = i;
                    res.push_back(pool->enqueue([this,di,party2]() {
                        io->recv_data(party2,di+party2,sizeof(int64_t));
                        io->flush(party2);
                    }));
                }
                joinNclean(res);
            }

            d = 0;
            for (int i = 1; i <= nP; i++) {

                d += di[i];
            #ifdef _debug_core_
                    cout<<"di:"<<di[i]<<endl;
            #endif
            }
            

            for (int i = 0; i < GC->cf->n1; i++)
            {
                delta[i] = ri_B[i].delta;

            }
            /*
            常量加法电路，<v>^B + d, 目前通过d被p1还原之后，比较容易做shareB，获得<d_i>^B，
            为了简化，在setup阶段，生成di_delta(依然是通过将除了p1之外的其他方设置为0，去简化通信量）
            */
            if (party == 1)
                prg.random_bool(delta+(GC->cf->n1),GC->cf->n2);
            else
                memset(delta+(GC->cf->n1),0,(GC->cf->n2)*sizeof(bool));

            memcpy(d_delta,delta+(GC->cf->n1),GC->cf->n2);

            GC->function_independent(delta);
	        // cout <<"FUNC_IND:\t"<<party<<"\n";	

            GC->function_dependent();
	        // cout <<"FUNC_DEP:\t"<<party<<"\n";

            delete[] delta;
            delta = nullptr;
            delete[] di;
            di = nullptr;
    }

    void A2B_setup(bool* res_delta, int64_t &d, bool* d_delta, ShareB ri_B[64], ShareA* r_A, ShareA *a, MPGC<nP>* GC) {

            ShareA ri[64];
            int64_t bais = 1;

            // int64_t di[nP+1] = {0};

            int num_in = (GC->cf->n1)+(GC->cf->n2);
            bool* delta = new bool[num_in];
            int64_t* di = new int64_t[nP+1];
            memset(di,0,sizeof(int64_t)*(nP+1));

            for (int i = 0; i < 64; i++) {
                ShareA tmp;
                GMW_A->randBit(ri+i);
                A2Bit(ri_B+i,ri+i);
                GMW_A->mulA_constant(&tmp,ri+i,bais);
                GMW_A->addA(r_A,r_A,&tmp);
                bais *= 2;
            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,ri+i);
                bool rec_B;
                GMW_B->rec_ShareB(rec_B,ri_B+i);
                cout<<(rec_A != (int64_t)rec_B ? "fail" : "success")<<endl;
            #endif
            }


            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,r_A);
                printf("rA:%lx\n",rec_A);
            #endif

            di[party] = a->delta - r_A->delta;

            if (party != 1) {
                io->send_data(1,di+party,sizeof(int64_t));
                io->flush(1);
            }
            else {
                vector<future<void>> res;//relic multi-thread problems...
                for (int i = 2; i <= nP; i++)
                {
                    int party2 = i;
                    res.push_back(pool->enqueue([this,di,party2]() {
                        io->recv_data(party2,di+party2,sizeof(int64_t));
                        io->flush(party2);
                    }));
                }
                joinNclean(res);
            }

            d = 0;
            for (int i = 1; i <= nP; i++) {

                d += di[i];
            #ifdef _debug_core_
                    cout<<"di:"<<di[i]<<endl;
            #endif
            }
            

            for (int i = 0; i < GC->cf->n1; i++)
            {
                delta[i] = ri_B[i].delta;

            }
            /*
            常量加法电路，<v>^B + d, 目前通过d被p1还原之后，比较容易做shareB，获得<d_i>^B，
            为了简化，在setup阶段，生成di_delta(依然是通过将除了p1之外的其他方设置为0，去简化通信量）
            */
            if (party == 1)
                prg.random_bool(delta+(GC->cf->n1),GC->cf->n2);
            else
                memset(delta+(GC->cf->n1),0,(GC->cf->n2)*sizeof(bool));

            memcpy(d_delta,delta+(GC->cf->n1),GC->cf->n2);

            GC->function_independent(delta);
	        // cout <<"FUNC_IND:\t"<<party<<"\n";	

            GC->function_dependent(res_delta);
	        // cout <<"FUNC_DEP:\t"<<party<<"\n";

            delete[] delta;
            delta = nullptr;
            delete[] di;
            di = nullptr;
    }
    
    void A2B_online(ShareB *res, int64_t d, bool* d_delta, ShareA* r_A, ShareB ri_B[64], ShareA *a, MPGC<nP>* GC) {

        int num_in = GC->cf->n1 + GC->cf->n2;
        int num_out = GC->cf->n3;

        int64_t D = 0;

        bool* input_value = new bool[num_in];

        for (int i = 0; i < GC->cf->n1; i++)
        {
            input_value[i] = ri_B[i].PublicVal;
            
        }

        D = a->PublicVal - r_A->PublicVal - d;
        for (int i = 0; i < GC->cf->n2; i++)
        {
            ShareB tmp;
            bool value = LSB(((uint64_t)D)>>i);
            GMW_B->set_ShareB(&tmp,value,1,d_delta+i);

            input_value[GC->cf->n1+i] = tmp.PublicVal;
        }

        bool* out_delta = new bool[num_out];
        bool* out_Public = new bool[num_out];
        memset(out_delta,0,sizeof(bool)*num_out);
        memset(out_Public,0,sizeof(bool)*num_out);

        GC->online(input_value,out_delta,out_Public);

        for (int i = 0; i < GC->cf->n3; i++)
        {
            res[i].delta = out_delta[i];
            res[i].PublicVal = out_Public[i];

            // cout<<"delta:"<<res[i].delta<<"\tpublicVal:"<<res[i].PublicVal<<endl;
        }

        delete[] input_value;
        delete[] out_delta;
        delete[] out_Public;
    }

    void A2B_online(ShareB *res, bool* out_delta, int64_t d, bool* d_delta, ShareA* r_A, ShareB ri_B[64], ShareA *a, MPGC<nP>* GC) {

        // time_t s,e;
        // s = clock();

        int num_in = GC->cf->n1 + GC->cf->n2;
        int num_out = GC->cf->n3;

        int64_t D = 0;

        bool* input_value = new bool[num_in];

        for (int i = 0; i < GC->cf->n1; i++)
        {
            input_value[i] = ri_B[i].PublicVal;
            
        }

        D = a->PublicVal - r_A->PublicVal - d;
        for (int i = 0; i < GC->cf->n2; i++)
        {
            ShareB tmp;
            bool value = LSB(((uint64_t)D)>>i);
            GMW_B->set_ShareB(&tmp,value,1,d_delta+i);

            input_value[GC->cf->n1+i] = tmp.PublicVal;
        }

        bool* out_Public = new bool[num_out];
        memset(out_Public,0,sizeof(bool)*num_out);
        // e = clock();
        // cout<<"prepare time: "<<((e-s)*1.0/CLOCKS_PER_SEC)<<"s"<<endl;

        // s = clock();
        GC->online(input_value,out_Public);
        // e = clock();
        // cout<<"GC_online time: "<<((e-s)*1.0/CLOCKS_PER_SEC)<<"s"<<endl;

        for (int i = 0; i < GC->cf->n3; i++)
        {
            res[i].delta = out_delta[i];
            res[i].PublicVal = out_Public[i];

            // cout<<"delta:"<<res[i].delta<<"\tpublicVal:"<<res[i].PublicVal<<endl;
        }

        delete[] input_value;
        delete[] out_Public;
    }
    
    void A2B_setup(bool* res_delta, bool* input_value, ShareA *a, MPGC<nP>* GC) {

            ShareA ri[64];
            ShareA r_A;
            ShareB ri_B[64];
            int64_t bais = 1;

            // int64_t di[nP+1] = {0};

            int num_in = (GC->cf->n1)+(GC->cf->n2);
            bool* delta = new bool[num_in];
            int64_t* di = new int64_t[nP+1];
            memset(di,0,sizeof(int64_t)*(nP+1));

            for (int i = 0; i < 64; i++) {
                ShareA tmp;
                GMW_A->randBit(ri+i);
                A2Bit(ri_B+i,ri+i);
                GMW_A->mulA_constant(&tmp,ri+i,bais);
                GMW_A->addA(&r_A,&r_A,&tmp);
                bais *= 2;
            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,ri+i);
                bool rec_B;
                GMW_B->rec_ShareB(rec_B,ri_B+i);
                cout<<(rec_A != (int64_t)rec_B ? "fail" : "success")<<endl;
            #endif
            }


            #ifdef _debug_core_
                int64_t rec_A;
                GMW_A->rec_ShareA(rec_A,r_A);
                printf("rA:%lx\n",rec_A);
            #endif

            di[party] = a->delta - r_A.delta;

            if (party != 1) {
                io->send_data(1,di+party,sizeof(int64_t));
                io->flush(1);
            }
            else {
                vector<future<void>> res;//relic multi-thread problems...
                for (int i = 2; i <= nP; i++)
                {
                    int party2 = i;
                    res.push_back(pool->enqueue([this,di,party2]() {
                        io->recv_data(party2,di+party2,sizeof(int64_t));
                        io->flush(party2);
                    }));
                }
                joinNclean(res);
            }

            int64_t d = 0;
            for (int i = 1; i <= nP; i++) {

                d += di[i];
            #ifdef _debug_core_
                    cout<<"di:"<<di[i]<<endl;
            #endif
            }
            

            for (int i = 0; i < GC->cf->n1; i++)
            {
                delta[i] = ri_B[i].delta;

            }

            /*
            常量加法电路，<v>^B + d, 目前通过d被p1还原之后，比较容易做shareB，获得<d_i>^B，
            为了简化，在setup阶段，生成di_delta(依然是通过将除了p1之外的其他方设置为0，去简化通信量）
            */
            if (party == 1)
                prg.random_bool(delta+(GC->cf->n1),GC->cf->n2);
            else
                memset(delta+(GC->cf->n1),0,(GC->cf->n2)*sizeof(bool));

            bool d_delta[64] = {0};
            memcpy(d_delta,delta+(GC->cf->n1),GC->cf->n2);

            int64_t D = 0;

            for (int i = 0; i < GC->cf->n1; i++)
            {
                input_value[i] = ri_B[i].PublicVal;
            
            }

            D = a->PublicVal - r_A.PublicVal - d;
            for (int i = 0; i < GC->cf->n2; i++)
            {
                ShareB tmp;
                bool value = LSB(((uint64_t)D)>>i);
                GMW_B->set_ShareB(&tmp,value,1,d_delta+i);

                input_value[GC->cf->n1+i] = tmp.PublicVal;
            }

            GC->function_independent(delta);
	        // cout <<"FUNC_IND:\t"<<party<<"\n";	

            GC->function_dependen_A2B(input_value,res_delta);
            // GC->function_dependent(res_delta);
	        // cout <<"FUNC_DEP:\t"<<party<<"\n";

            delete[] delta;
            delta = nullptr;
            delete[] di;
            di = nullptr;
    }

    void A2B_online(ShareB *res, bool* out_delta, bool input_value[64],  ShareA *a, MPGC<nP>* GC) {

        // time_t s,e;
        struct timeval t1,t2;
        double timeuse;
        
        // s = clock();

        
        int num_out = GC->cf->n3;
        bool* out_Public = new bool[num_out];
        memset(out_Public,0,sizeof(bool)*num_out);
        // e = clock();
        // cout<<"prepare time: "<<((e-s)*1.0/CLOCKS_PER_SEC)<<"s"<<endl;

        // s = clock();
        gettimeofday(&t1,NULL);
        GC->online_A2B(input_value,out_Public);
        gettimeofday(&t2,NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        // e = clock();
        cout<<"GC_online throughput: "<<1.0/timeuse<<"opt/s"<<endl;

        for (int i = 0; i < GC->cf->n3; i++)
        {
            res[i].delta = out_delta[i];
            res[i].PublicVal = out_Public[i];

            // cout<<"delta:"<<res[i].delta<<"\tpublicVal:"<<res[i].PublicVal<<endl;
        }

        delete[] out_Public;
    }

    void A2B_setup_api(A2B_ctx2* ctx, ShareA *a, MPGC<nP>* &GC) {
        
        // GC = new MPGC<nP>(ios, pool, party, &cf);
        A2B_setup(ctx->res_delta, ctx->input_value, a, GC);

    }
    void A2B_online_api(ShareB *res, A2B_ctx2* ctx, ShareA *a, MPGC<nP>* &GC) {

        A2B_online(res, ctx->res_delta, ctx->input_value, a, GC);
        // delete GC;
    }

    void A2B_setup_api(A2B_ctx* ctx, ShareA *a, MPGC<nP>* &GC) {
        
        // GC = new MPGC<nP>(ios, pool, party, &cf);
        A2B_setup(ctx->res_delta, ctx->d, ctx->d_delta, ctx->ri_B, &(ctx->r_A), a, GC);

    }
    void A2B_online_api(ShareB *res, A2B_ctx* ctx, ShareA *a, MPGC<nP>* &GC) {

        A2B_online(res, ctx->res_delta, ctx->d, ctx->d_delta, &(ctx->r_A), ctx->ri_B, a, GC);
        // delete GC;
    }


    void scalarP_setup(int64_t &delta_c, int64_t* delta_ab, int64_t* a_delta, int64_t* b_delta, int length) {

        prg.random_data(&delta_c,sizeof(int64_t)); // generate delta_c

        // #pragma omp parallel for //并行的时候会出错，目前还不知道问题在哪里，后续可以进行优化
        for (int i = 0; i < length; i++)
        {
            GMW_A->setupMULT(delta_ab[i],a_delta[i],b_delta[i]);
        }
    }

    void scalarP_online(ShareA* res, int64_t &delta_c, int64_t* delta_ab, ShareA* A, ShareA* B, int length) {

        int64_t Delta_share = 0;
        int64_t Delta_sum = 0;
        if (party == 1) {
            for (int i = 0; i < length; i++)
            {
                Delta_share = Delta_share + (A[i].PublicVal * B[i].PublicVal - A[i].PublicVal * B[i].delta - B[i].PublicVal * A[i].delta + delta_ab[i]);
            }

            Delta_share = Delta_share + delta_c;
        }
        else {
            for (int i = 0; i < length; i++)
            {
                Delta_share = Delta_share + delta_ab[i] - A[i].PublicVal * B[i].delta - B[i].PublicVal * A[i].delta;
            }

            Delta_share = Delta_share + delta_c;
        }

        GMW_A->open(Delta_sum,Delta_share);

        res->delta = delta_c;
        res->PublicVal = Delta_sum;  
    }

    //more slowly
    void scalarP_online_digen(ShareA* res, int64_t &delta_c, int64_t* delta_ab, int64_t* Delta_a, int64_t* delta_a, int64_t* Delta_b, int64_t* delta_b, int length) {

        int64_t Delta_share = 0;
        int64_t Delta_sum = 0;

        // RowVectorX<int64_t> PublicVal_A(length) , PublicVal_B(length), Delta_A(length) , Delta_B(length), Delta_ab(length);

        time_t s,e;


        Map<RowVectorX<int64_t>> PublicVal_A(Delta_a,length);
        Map<RowVectorX<int64_t>> PublicVal_B(Delta_b,length);
        Map<RowVectorX<int64_t>> Delta_A(delta_a,length);
        Map<RowVectorX<int64_t>> Delta_B(delta_b,length);
        Map<RowVectorX<int64_t>> Delta_ab(delta_ab,length);
    
        // for (int i = 0; i < length; i++)
        // {
        //     PublicVal_A.col(i) << A[i].PublicVal;
        //     PublicVal_B.col(i)  << B[i].PublicVal;
        //     Delta_A.col(i)  << A[i].delta;
        //     Delta_B.col(i)  << B[i].delta;
        //     Delta_ab.col(i)  << delta_ab[i];
        // }
        
        s = clock();
        if (party == 1)
            Delta_share = (PublicVal_A * PublicVal_B.transpose() - PublicVal_A * Delta_B.transpose() - PublicVal_B * Delta_A.transpose() + RowVectorX<int64_t>::Ones(length) * Delta_ab.transpose()).value() + delta_c;
        else
            Delta_share = (RowVectorX<int64_t>::Ones(length) * Delta_ab.transpose() - PublicVal_A * Delta_B.transpose() - PublicVal_B * Delta_A.transpose()).value() + delta_c;
        
        GMW_A->open(Delta_sum,Delta_share);

        e = clock();
        cout<<"interl time: "<<((double)(e-s)*1.0/CLOCKS_PER_SEC)<<"s"<<endl;

        res->delta = delta_c;
        res->PublicVal = Delta_sum;  
    }

    void MM_setup(Eigen::MatrixX<int64_t> &delta_C, Eigen::MatrixX<int64_t> &delta_AB, Eigen::MatrixX<int64_t> delta_a, Eigen::MatrixX<int64_t> delta_b) {

        if(delta_a.cols()!=delta_b.rows()) return;


        delta_C.setRandom();
        RowVectorX<int64_t> temp_ab(delta_a.cols());

        int64_t aa = 0;
        int64_t bb = 0;
        int64_t cc = 0;

        for (int i = 0; i < delta_a.rows(); i++)
        {
            for (int j = 0; j < delta_b.cols(); j++)
            {
                temp_ab.setZero();
                for (int k = 0; k < delta_a.cols(); k++)
                {   aa = delta_a(i,k);
                    bb = delta_b(k,j);
                    GMW_A->setupMULT(cc,aa,bb);
                    temp_ab.col(k) << cc;
                }
                delta_AB(i,j) = (RowVectorX<int64_t>::Ones(delta_a.cols()) * temp_ab.transpose());
            }
        }
    

    }

    void MM_online(Eigen::MatrixX<int64_t> &PublicVal_C, Eigen::MatrixX<int64_t> delta_C, Eigen::MatrixX<int64_t> &delta_AB, 
    Eigen::MatrixX<int64_t> PublicVal_a, Eigen::MatrixX<int64_t> PublicVal_b, Eigen::MatrixX<int64_t> delta_a, Eigen::MatrixX<int64_t> delta_b) {

        Eigen::MatrixX<int64_t> temp_C(PublicVal_C.rows(),PublicVal_C.cols());
        temp_C.setZero();
        int64_t Delta_C[PublicVal_C.size()] = {0};
        
        if (party == 1)
            temp_C = PublicVal_a*PublicVal_b - PublicVal_a*delta_b - delta_a*PublicVal_b + delta_AB + delta_C;
        else
            temp_C = delta_AB + delta_C - PublicVal_a*delta_b - delta_a*PublicVal_b;


        
        GMW_A->open_vec(Delta_C,temp_C.data(),temp_C.size());

        Map<MatrixX<int64_t>> temp_C1(Delta_C,PublicVal_C.rows(),PublicVal_C.cols());

        PublicVal_C = temp_C1;
    }


    void MMTR_setup(Eigen::MatrixX<int64_t> &Public_R_tr, Eigen::MatrixX<int64_t> &delta_R_tr, 
    Eigen::MatrixX<int64_t> &delta_R, Eigen::MatrixX<int64_t> &delta_AB, 
    Eigen::MatrixX<int64_t> delta_a, Eigen::MatrixX<int64_t> delta_b, int truncation) {

        if(delta_a.cols()!=delta_b.rows()) 
        {   cout<<"error"<<endl;
            return;
        }

        MatrixX<int64_t> Public_R(delta_R.rows(),delta_R.cols());
        Public_R.setZero();
        delta_R.setZero();
        delta_R_tr.setZero();
        Public_R_tr.setZero();
        int64_t bais = 1;
        int64_t bais2 = 1;
        

        for (int i = 0; i < 64; i++) {
                int64_t Public_tmp_val[delta_R.size()] = {0};
                int64_t Delta_tmp_val[delta_R.size()] = {0};
                GMW_A->randBit(Public_tmp_val,Delta_tmp_val,delta_R.size());
                Map<MatrixX<int64_t>> Public_tmp(Public_tmp_val,delta_R.rows(),delta_R.cols());
                Map<MatrixX<int64_t>> Delta_tmp(Delta_tmp_val,delta_R.rows(),delta_R.cols());

                //Combine into R
                MatrixX<int64_t> tmp(delta_R.rows(),delta_R.cols());
                tmp.setZero();
                GMW_A->mulA_constant_setup_mat(tmp,Delta_tmp,bais);
                GMW_A->addA_setup_mat(delta_R,delta_R,tmp);
                tmp.setZero();
                GMW_A->mulA_constant_online_mat(tmp,Public_tmp,bais);
                GMW_A->addA_online_mat(Public_R,Public_R,tmp);
                bais *= 2;
                #ifdef _debug_core_

                    int64_t tmp_value[delta_R.size()] = {0};
                    GMW_A->rec_ShareA_vec(tmp_value,Public_tmp_val,Delta_tmp_val,delta_R.size());
                    int64_t R_value[delta_R.size()] = {0};
                    GMW_A->rec_ShareA_vec(R_value,Public_R.data(),delta_R.data(),Public_R.size());
                    cout<<"R:";
                    for (int i = 0; i < Public_R.size(); i++)
                    {
                        cout<<R_value[i]<<std::hex<<"   ";
                        // cout<<(R_value[i]== tmp_value[i] ?"success":"fail")<<" ";
                    }
                    cout<<tmp_value[Public_R.size()-1]<<endl;
                #endif
                //Combine into R_trunc
                if(i >= truncation) {
                    MatrixX<int64_t> tmp2(delta_R_tr.rows(),delta_R_tr.cols());
                    tmp2.setZero();
                    GMW_A->mulA_constant_setup_mat(tmp2,Delta_tmp,bais2);
                    GMW_A->addA_setup_mat(delta_R_tr,delta_R_tr,tmp2);
                    tmp2.setZero();
                    GMW_A->mulA_constant_online_mat(tmp2,Public_tmp,bais2);
                    GMW_A->addA_online_mat(Public_R_tr,Public_R_tr,tmp2);
                    bais2 *= 2;
                    #ifdef _debug_core_
                    int64_t Rtr_value[delta_R_tr.size()] = {0};
                    GMW_A->rec_ShareA_vec(Rtr_value,Public_R_tr.data(),delta_R_tr.data(),delta_R_tr.size());
                    cout<<"R_trunc:";
                    for (int i = 0; i < delta_R_tr.size(); i++)
                    {
                        cout<<Rtr_value[i]<<std::hex<<"   ";
                        // cout<<(R_value[i]== tmp_value[i] ?"success":"fail")<<" ";
                    }
                    cout<<tmp_value[Public_R.size()-1]<<endl;
                    #endif
                }
                if(i == 63) {
                    for (int j = 0; j < truncation; j++)
                    {
                        MatrixX<int64_t> tmp2(delta_R_tr.rows(),delta_R_tr.cols());
                        tmp2.setZero();
                        GMW_A->mulA_constant_setup_mat(tmp2,Delta_tmp,bais2);
                        GMW_A->addA_setup_mat(delta_R_tr,delta_R_tr,tmp2);
                        tmp2.setZero();
                        GMW_A->mulA_constant_online_mat(tmp2,Public_tmp,bais2);
                        GMW_A->addA_online_mat(Public_R_tr,Public_R_tr,tmp2);
                        bais2 *= 2;
                    }
                }
            }
            #ifdef _debug_core_
                    int64_t Rtr_value[delta_R_tr.size()] = {0};
                    GMW_A->rec_ShareA_vec(Rtr_value,Public_R_tr.data(),delta_R_tr.data(),delta_R_tr.size());
                    int64_t R_value[delta_R.size()] = {0};
                    GMW_A->rec_ShareA_vec(R_value,Public_R.data(),delta_R.data(),Public_R.size());
                    cout<<"R:";
                    for (int i = 0; i < delta_R_tr.size(); i++)
                    {
                        // cout<<Rtr_value[i]<<std::hex<<endl;
                        // cout<<R_value[i]<<std::hex<<endl;
                        cout<<(Rtr_value[i] == (R_value[i]>>truncation) ? "success" : "fail")<<endl;
                        // cout<<(R_value[i]== tmp_value[i] ?"success":"fail")<<" ";
                    }
            #endif
            
        //generate delta_AB
        RowVectorX<int64_t> temp_ab(delta_a.cols());

        int64_t aa = 0;
        int64_t bb = 0;
        int64_t cc = 0;

        for (int i = 0; i < delta_a.rows(); i++)
        {
            for (int j = 0; j < delta_b.cols(); j++)
            {
                temp_ab.setZero();
                for (int k = 0; k < delta_a.cols(); k++)
                {   aa = delta_a(i,k);
                    bb = delta_b(k,j);
                    GMW_A->setupMULT(cc,aa,bb);
                    temp_ab.col(k) << cc;
                }
                delta_AB(i,j) = (RowVectorX<int64_t>::Ones(delta_a.cols()) * temp_ab.transpose());
            }
        }

        //

        if(party == 1)
            delta_R = Public_R - delta_R;
        else
            delta_R = -delta_R;
    }

    void MMTR_online(Eigen::MatrixX<int64_t> &Public_C, Eigen::MatrixX<int64_t> &delta_C, 
    Eigen::MatrixX<int64_t> Public_R_tr, Eigen::MatrixX<int64_t> delta_R_tr, 
    Eigen::MatrixX<int64_t> delta_R, Eigen::MatrixX<int64_t> delta_AB, 
    Eigen::MatrixX<int64_t> Puiblic_a, Eigen::MatrixX<int64_t> delta_a, Eigen::MatrixX<int64_t> Public_b, Eigen::MatrixX<int64_t> delta_b,  
    int truncation)  {
        
        int64_t Z_val[delta_R.size()] = {0};
        MatrixX<int64_t> Z(delta_R.rows(),delta_R.cols());
        Z.setZero();

        if (party == 1)
            Z = Puiblic_a*Public_b - Puiblic_a*delta_b - delta_a*Public_b + delta_AB;
        else
            Z = delta_AB - Puiblic_a*delta_b - delta_a*Public_b;

        Z = Z + delta_R;

        GMW_A->open_vec(Z_val,Z.data(),Z.size());

        //turncate Z locally
        Map<MatrixX<int64_t>> Z_mat(Z_val,delta_R.rows(),delta_R.cols());

        Public_C = (Z_mat / (1 << truncation)) - Public_R_tr;
        delta_C = -delta_R_tr;
    }


    ~CORE_protocol() {
        delete GMW_A;
        delete GMW_B;

    }
};

#endif