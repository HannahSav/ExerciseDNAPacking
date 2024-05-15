#pragma once
#include <vector>
#include <string>


struct float3{
	float x,y,z;

    float length() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    float3 distance(float3 other) {
        return {x - other.x, y - other.y, z - other.z};
    }

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
};

void saveCoordinates(const std::string filename,
                     const std::string modifier,
                    //  const std::vector<Atom>& atoms,
                     const std::vector<float3>& r,
                     const std::vector<float3>& v, 
                     int& N);



float3 f_i_new(int& i, int& N_pairs, HiCData& hic, int& max_U_ij, std::vector<float3>& r);

float3 v_i_new(float3 v_prev, float3 f);

float3 r_i_new(float3 r_prev, float3 v_next);

std::pair<std::vector<float3>, std::vector<float3> > step(int& N, std::vector<float3>& v, std::vector<float3>& r, std::vector<std::vector<int> >& U, int& max_U_ij );

std::vector<float3> counting(HiCData& hic);


