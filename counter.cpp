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
#define NSTEPS 10 //1000000
#define STRIDE 1
#define L 1.0


void saveCoordinates(const std::string filename,
                     const std::string modifier,
                     const std::vector<float3>& r,
                     const std::vector<float3>& v, 
                     int& N)
{
    FILE* fout = fopen(filename.c_str(), modifier.c_str());
    fprintf(fout, "DNA\n%d\n", N);
    for (int i = 0; i < N; i++)
    {
        fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
            i, "DNA",
            "K", i, 
            r[i].x, r[i].y, r[i].z,
            v[i].x, v[i].y, v[i].z
        );
    }
    fprintf(fout, "%8.3f %8.3f %8.3f\n", L, L, L);
    fclose(fout);
}

float3 pair_grad(HarmonicPair& pair, float3& r_i, float3& r_j){
    float3 r_ij = r_i - r_j;
    float3 grad_res = -(r_ij*2*pair.kc - r_ij/r_ij.length()*pair.kc*pair.r0);
    return grad_res;
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

std::pair<std::vector<float3>, std::vector<float3> > step(std::vector<HarmonicPair> all_pairs, HiCData& hic, std::vector<float3>& v, std::vector<float3>& r, int& max_U_ij ){
    std::vector<float3> r_new(hic.atomCount, {0.0f, 0.0f, 0.0f});
    std::vector<float3> v_new(hic.atomCount, {0.0f, 0.0f, 0.0f});
    std::vector<float3> f_new(hic.atomCount, {0.0f, 0.0f, 0.0f});

    int N = hic.atomCount;
    int N_pairs = hic.pairCount;

    float3 fi;
    float3 fj;
    HarmonicPair pair;
    for(int p = 0; p < all_pairs.size(); p++){
        pair = all_pairs[p];
        printf("\n%d - pair -%d\n", pair.i, pair.j);
        fi = pair_grad(pair, r[pair.i], r[pair.j]); //градиент пары
        printf("\n Grad res %f %f %f \n", fi.x, fi.y, fi.z);
        fj = -fi;
        f_new[pair.i] += fi;
        f_new[pair.j] += fj;
    }

    printf("\nCounted f");
    float3 v_res;
    for (int i = 0; i < N; i++){
        v_res = v_i_new(v[i], f_new[i]);
        v_new.push_back(v_res);
        printf("\n V res %f %f %f \n", v_res.x, v_res.y, v_res.z);
    }
    printf("\nCounted v");
    float3 r_res;
    for (int i = 0; i < N; i++){
        r_res = r_i_new(r[i], v_new[i]);
        r_new.push_back(r_res);
        printf("\n R res %f %f %f \n", r_res.x, r_res.y, r_res.z);
    }
    printf("\nCounted r\n");
    return std::make_pair(v_new, r_new);
}

std::vector<float3> counting(HiCData& hic){
    //counting start meanings
    std::vector<float3> r_start;
    std::vector<float3> v_start(hic.atomCount, {0.0f, 0.0f, 0.0f});

    for(int i = 0; i < hic.atomCount; i++){
        r_start.push_back(float3{1.0f*i - hic.atomCount/2, 0.0, 0.0});
    }

    int max_U_ij = 0;
    for (int k = 0; k < hic.pairCount; ++k){
        if (hic.pairs[k].n > max_U_ij){
            max_U_ij = hic.pairs[k].n;
        }
    }
    

    //harmonicPair (AllPairs)
    std::vector<HarmonicPair> all_pairs;
    HarmonicPair one_pair;
    //adiing V_cov
    for (int i = 0; i < hic.atomCount-1; i++){
        one_pair.i = i;
        one_pair.j = i+1;
        one_pair.kc = k_c;
        one_pair.r0 = r_0;
        all_pairs.push_back(one_pair);
    }

    //adding V_hic
    for (int p = 0; p < hic.pairCount; p++){
        one_pair.i = hic.pairs[p].i;
        one_pair.j = hic.pairs[p].j;
        one_pair.kc = 50.0 * hic.pairs[p].n/max_U_ij;
        one_pair.r0 = 0.19 + 0.01*max_U_ij/hic.pairs[p].n;
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
        step_res = step(all_pairs, hic, v_new, r_new, max_U_ij);
        v_new = step_res.first;
        r_new = step_res.second;
        // for(int i = 0; i < hic.atomCount; i++){
        //     temp += m*(v_new[i].x*v_new[i].x + v_new[i].y*v_new[i].y + v_new[i].z*v_new[i].z);
        // }
    }
    return r_new;
} 