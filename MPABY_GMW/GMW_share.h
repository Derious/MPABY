// #ifndef GMWSHAREA_H__
// #define GMWSHAREA_H__

// #include "../MPABY_GC/GC_netmp.h"
// #include "../MPABY_GC/GC_helper.h"
// using namespace emp;

// #define _debug
// template<int nP>
// class GMWShareA
// {
// public:

//     NetIOMP<nP> *io;
// 	ThreadPool * pool;
// 	int party;
//     PRG prg;
//     int64_t PublicValue;
//     int64_t delta;

//     GMWShareA(NetIOMP<nP>* io, ThreadPool * pool, int party) {
//         this->io = io;
//         this->pool = pool;
//         this->party = party;
//     }

//     void ShareA(int64_t value, int p) {

//             vector<future<void>> res;//relic multi-thread problems...

//             if(p == party) {
//                 prg.random_data(&delta,8);
//                 PublicValue = delta + value;
//                 for (int i = 1; i <= nP; i++) {
//                     int party2 = i;
//                     if (i != party) {
//                         res.push_back(pool->enqueue([this,party2]() {
//                             io->send_data(party2,&PublicValue,sizeof(int64_t));
//                             io->flush(party2);
// 				            }));
//                     }
//                 }
//                 joinNclean(res);
//             }
//             else {
//                 delta = 0;
//                 io->recv_data(p,&PublicValue,sizeof(int64_t));
//                 io->flush(p);
//             }

// #ifdef _debug

//             cout<<"Public Delta:"<<PublicValue<<endl;
//             cout<<"Share delta:"<<delta<<endl;

// #endif
                
//     }

//     void RecA(int64_t &value){
        
//         int64_t tmp = 0;
//         int64_t *delta_i = new int64_t[nP+1];
//         delta_i[party] = delta;

//         vector<future<void>> res;//relic multi-thread problems...
//         for (int i = 1; i <= nP; i++)
//         {   
//             int party2 = i;
//             if (i!=party)
//             {
//                 res.push_back(pool->enqueue([this,party2]() {
//                     io->send_data(party2,&delta,sizeof(int64_t));
//                     io->flush(party2);
// 				}));
//                 res.push_back(pool->enqueue([this,delta_i,party2]() {
//                     io->recv_data(party2,delta_i+party2,sizeof(int64_t));
//                     io->flush(party2);
//                 }));
//             }
//         }
//         joinNclean(res);

//         for (int i = 1; i <= nP; i++)
//         {
//             tmp += delta_i[i];
//         }
//         value = PublicValue - tmp;
        
// #ifdef _debug
//         cout<<"value:"<<value<<endl;
// #endif
//     }


//     ~GMWShareA() {

//     }
// };



// template<int nP>
// class GMWShareB
// {
// public:

//     NetIOMP<nP> *io;
// 	ThreadPool * pool;
// 	int party;
//     PRG prg;
//     bool PublicValue;
//     bool delta;

//     GMWShareB(NetIOMP<nP>* io, ThreadPool * pool, int party) {
//         this->io = io;
//         this->pool = pool;
//         this->party = party;
//     }

//     void ShareB(bool value, int p) {

//             vector<future<void>> res;//relic multi-thread problems...

//             if(p == party) {
//                 prg.random_bool(&delta,1);
//                 PublicValue = delta != value;
//                 for (int i = 1; i <= nP; i++) {
//                     int party2 = i;
//                     if (i != party) {
//                         res.push_back(pool->enqueue([this,party2]() {
//                             io->send_data(party2,&PublicValue,sizeof(bool));
//                             io->flush(party2);
// 				            }));
//                     }
//                 }
//                 joinNclean(res);
//             }
//             else {
//                 delta = 0;
//                 io->recv_data(p,&PublicValue,sizeof(bool));
//                 io->flush(p);
//             }

// #ifdef _debug

//             cout<<"Public Delta:"<<PublicValue<<endl;
//             cout<<"Share delta:"<<delta<<endl;

// #endif
                
//     }
//     void RecB(bool &value){
        
//         bool tmp = 0;
//         bool *delta_i = new bool[nP+1];
//         delta_i[party] = delta;

//         vector<future<void>> res;//relic multi-thread problems...
//         for (int i = 1; i <= nP; i++)
//         {   
//             int party2 = i;
//             if (i!=party)
//             {
//                 res.push_back(pool->enqueue([this,party2]() {
//                     io->send_data(party2,&delta,sizeof(bool));
//                     io->flush(party2);
// 				}));
//                 res.push_back(pool->enqueue([this,delta_i,party2]() {
//                     io->recv_data(party2,delta_i+party2,sizeof(bool));
//                     io->flush(party2);
//                 }));
//             }
//         }
//         joinNclean(res);

//         for (int i = 1; i <= nP; i++)
//         {
//             tmp = tmp != delta_i[i];
//         }
//         value = PublicValue != tmp;
        
// #ifdef _debug
//         cout<<"value:"<<value<<endl;
// #endif
//     }

//     ~GMWShareB() {

//     }
// };

// #endif GMWSHAREA_H__