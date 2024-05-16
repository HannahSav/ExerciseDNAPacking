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
#define T 50 //K
#define NSTEPS 10000 //1000000
#define STRIDE 100
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
    float3 r_ij = - r_i + r_j;
    float3 grad_res = -(r_ij*2*pair.kc - r_ij/r_ij.length()*pair.kc*2*pair.r0);
    return grad_res;
}

float3 v_i_new(float3 v_prev, float3 f, std::mt19937& gen, std::normal_distribution<>& d){
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

void step(const std::vector<HarmonicPair>& all_pairs, const HiCData& hic, std::vector<float3>& v, std::vector<float3>& r, std::vector<float3>& f, int max_U_ij, std::mt19937& gen, std::normal_distribution<>& d){
    int N = hic.atomCount;

    float3 fi;
    float3 fj;
    HarmonicPair pair;

    // Обнуление вектора сил
    for (int i = 0; i < N; ++i)
        f[i] = {0.0f, 0.0f, 0.0f};

    for(int p = 0; p < all_pairs.size(); p++){
        pair = all_pairs[p];
        fi = pair_grad(pair, r[pair.i], r[pair.j]); 
        fj = -fi;
        f[pair.i] += fi;
        f[pair.j] += fj;
    }

    float3 v_res;
    for (int i = 0; i < N; i++){
        v_res = v_i_new(v[i], f[i], gen, d);
        v[i] = v_res;
    }

    float3 r_res;
    for (int i = 0; i < N; i++){
        r_res = r_i_new(r[i], v[i]);
        r[i] = r_res;
    }
}

std::vector<float3> counting(HiCData& hic){
    std::vector<float3> r_start;
    std::vector<float3> v_start(hic.atomCount, {0.0f, 0.0f, 0.0f});
    std::vector<float3> f_start(hic.atomCount, {0.0f, 0.0f, 0.0f});

    for(int i = 0; i < hic.atomCount; i++){
        r_start.push_back(float3{1.0f*i - hic.atomCount/2, 0.0, 0.0});
    }

    int max_U_ij = 0;
    for (int k = 0; k < hic.pairCount; ++k){
        if (hic.pairs[k].n > max_U_ij){
            max_U_ij = hic.pairs[k].n;
        }
    }
    
    std::vector<HarmonicPair> all_pairs;
    HarmonicPair one_pair;

    for (int i = 0; i < hic.atomCount-1; i++){
        one_pair.i = i;
        one_pair.j = i+1;
        one_pair.kc = k_c* 10.0;
        one_pair.r0 = r_0;
        all_pairs.push_back(one_pair);
    }

    for (int p = 0; p < hic.pairCount; p++){
        one_pair.i = hic.pairs[p].i;
        one_pair.j = hic.pairs[p].j;
        one_pair.kc = 50.0 * hic.pairs[p].n/max_U_ij * 10.0;
        one_pair.r0 = 0.19 + 0.01*max_U_ij/hic.pairs[p].n;
        all_pairs.push_back(one_pair);
    }

    std::vector<float3> r_new = r_start;
    std::vector<float3> v_new = v_start;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);

    float temp = 0.0;
    std::string mod = "w";
    for(int t = 0; t < NSTEPS; t++){
        if (t % STRIDE == 0){
            printf("%d %f\n", t, temp);
            saveCoordinates("coord.gro", mod, r_new, v_new, hic.atomCount);
            mod = "a";
        }
        temp = 0.0;
        step(all_pairs, hic, v_new, r_new, f_start, max_U_ij, gen, d);
    }
}
