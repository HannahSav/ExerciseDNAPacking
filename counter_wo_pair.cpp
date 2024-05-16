#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <random>
#include "hicio.h"
#include "counter.h"

#define m 1.0
#define r_0 0.1
#define k_c 100
#define lambda 50
#define k_b 8.31 * 0.001 //кДж/К*моль
#define tau 0.001 // пикасекунды
#define T 300 //K
#define NSTEPS 500 //1000000
#define STRIDE 10
#define L 1.0


void saveCoordinates(const std::string filename,
                     const std::string modifier,
                    //  const std::vector<Atom>& atoms,
                     const std::vector<float3>& r,
                     const std::vector<float3>& v, 
                     int& N)
{
    FILE* fout = fopen(filename.c_str(), modifier.c_str());
    fprintf(fout, "DNA\n%d\n", N);
    for (int i = 0; i < N; i++)
    {
        // printf(" mem %f", r[i].x);
        fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
            i, "DNA",
            "K", i, 
            // atoms[i].name.c_str(), atoms[i].index,
            r[i].x, r[i].y, r[i].z,
            v[i].x, v[i].y, v[i].z
        );
    }
    fprintf(fout, "%8.3f %8.3f %8.3f\n", L, L, L);
    fclose(fout);
}


float3 f_i_new(int& i, int& N_pairs, std::vector<std::vector<int> >& k_for_i, HiCData& hic, int& max_U_ij, std::vector<float3>& r){
    float3 V_cov = {0, 0, 0};
    float3 V_hic = {0, 0, 0};
    
    if (i != 0){
        float3 dist1 = r[i].distance(r[i-1]);
        V_cov -= dist1 * k_c*2 - dist1/ dist1.length() * 2 * k_c *r_0;
    }
    if (i != hic.atomCount - 1){
        float3 dist2 = r[i+1].distance(r[i]);
        V_cov -= dist2 * k_c*2 - dist2/dist2.length() * 2 * k_c *r_0;
    }
    
    
    
    float k_ij;
    float r_ij_0;

     for(int k = 0; k < k_for_i[i].size(); k++){
        int num_of_pair = k_for_i[i][k];
        int i_pair = hic.pairs[num_of_pair].i;
        int j_pair = hic.pairs[num_of_pair].j;
        int n_pair = hic.pairs[num_of_pair].n;
        float3 dist0 = r[i_pair].distance(r[j_pair]);
        k_ij = 50*n_pair/max_U_ij;
        r_ij_0 = 0.19 + 0.01*max_U_ij/n_pair;
        V_hic -= dist0*k_ij * 2 - dist0/dist0.length()*2*r_ij_0;
    }
    return V_cov+V_hic;
}

float3 v_i_new(float3 v_prev, float3 f){
    //Creating normal distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);
    float3 r_i_f;

    r_i_f.x = d(gen);
    r_i_f.y = d(gen);
    r_i_f.z = d(gen);
    
    float3 n_comp;
    float3 v;

    n_comp = r_i_f*std::sqrt(2.0*k_b*T*lambda*m/tau);
    v = (f + v_prev*m/tau - v_prev*lambda*m/2 + n_comp)/(m/tau + lambda*m/2);

    return v;
}

float3 r_i_new(float3 r_prev, float3 v_next){
    return v_next*tau + r_prev;
}

std::pair<std::vector<float3>, std::vector<float3> > step(HiCData& hic, std::vector<std::vector<int> >& k_for_i, std::vector<float3>& v, std::vector<float3>& r, int& max_U_ij ){
    std::vector<float3> r_new;
    std::vector<float3> v_new;
    std::vector<float3> f_new;

    int N = hic.atomCount;
    int N_pairs = hic.pairCount;
    for(int i = 0; i < N; i++){
        float3 res = f_i_new(i, N_pairs, k_for_i, hic, max_U_ij, r);
        f_new.push_back(res);
    }
    printf("\nCounted f");
    for (int i = 0; i < N; i++){
        v_new.push_back(v_i_new(v[i], f_new[i]));
    }
    printf("\nCounted v");
    for (int i = 0; i < N; i++){
        r_new.push_back(r_i_new(r[i], v_new[i]));
    }
    printf("\nCounted r\n");
    return std::make_pair(v_new, r_new);
}

std::vector<float3> counting(HiCData& hic){
    //counting start meanings
    std::vector<float3> r_start;
    std::vector<float3> v_start(hic.atomCount, {0.0f, 0.0f, 0.0f});

    for(int i = 0; i < hic.atomCount; i++){
        r_start.push_back(float3{1.0f*i, 0.0, 0.0});
        // printf("kek %f\n", 0.1*i);
    }

    int max_U_ij = 0;
    for (int k = 0; k < hic.pairCount; ++k){
        if (hic.pairs[k].n > max_U_ij){
            max_U_ij = hic.pairs[k].n;
        }
    }
    


    //numbers of pairs for i
    std::vector<std::vector<int> > k_for_i;
    std::vector<int> one_k_for_i;
    for (int i = 0; i < hic.atomCount; i++){
        one_k_for_i.clear();
        for(int k = 0; k < hic.pairCount; k++){
            int i_pair = hic.pairs[k].i;
            int j_pair = hic.pairs[k].j;
            int n_pair = hic.pairs[k].n;
            if (i_pair == i || j_pair == i){
                one_k_for_i.push_back(k);
            }
        }
        k_for_i.push_back(one_k_for_i);
    }


    //iterations by time
    std::vector<float3> r_new = r_start;
    std::vector<float3> v_new = v_start;
    std::pair<std::vector<float3>, std::vector<float3> > step_res;

    float temp = 0.0;
    std::string mod = "w";
    for(int t = 0; t < NSTEPS; t++){
        if (t % STRIDE == 0){
            printf("%d %f\n", t, temp);
            saveCoordinates("coord.gro", mod, r_new, v_new, hic.atomCount);
            mod = "a";
        }
        temp = 0.0;
        step_res = step(hic, k_for_i, v_new, r_new, max_U_ij);
        v_new = step_res.first;
        r_new = step_res.second;
        // for(int i = 0; i < hic.atomCount; i++){
        //     temp += m*(v_new[i].x*v_new[i].x + v_new[i].y*v_new[i].y + v_new[i].z*v_new[i].z);
        // }
    }
    return r_new;
} 