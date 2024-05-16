#pragma once
#include <vector>
#include <string>


struct HarmonicPair{
    int i, j;
    float kc, r0;
};

struct float3{
	float x,y,z;

    float length() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    // float3 distance(float3 other) {
    //     return {x - other.x, y - other.y, z - other.z};
    // }

    float3 operator/(float scalar) const{
        return {x/scalar, y/scalar, z/scalar};
    }
    float3 operator*(float scalar) const{
        return {x * scalar, y*scalar, z*scalar};
    }
    float3 operator+(float3 other) const{
        return{x+other.x, y+other.y, z+other.z};
    }
    float3 operator+=(float3 other) const{
        return{x+other.x, y+other.y, z+other.z};
    }
    float3 operator-(float3 other) const{
        return{x-other.x, y-other.y, z-other.z};
    }
    float3 operator-() const {
        return {-x, -y, -z};
    }
};

void saveCoordinates(const std::string filename,
                     const std::string modifier,
                    //  const std::vector<Atom>& atoms,
                     const std::vector<float3>& r,
                     const std::vector<float3>& v, 
                     int& N);


float3 pair_grad(HarmonicPair& pair, float3& r_i, float3& r_j);

float3 v_i_new(float3 v_prev, float3 f);

float3 r_i_new(float3 r_prev, float3 v_next);

std::pair<std::vector<float3>, std::vector<float3> > step(std::vector<HarmonicPair>& all_pairs, HiCData& hic, std::vector<float3>& v, std::vector<float3>& r, int& max_U_ij );

std::vector<float3> counting(HiCData& hic);