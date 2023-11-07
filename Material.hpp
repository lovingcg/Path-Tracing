//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE, MICROFACET, MIRROR};

class Material{
private:
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    float DistributionGGX(float NdotH, float roughness) {
        float a = roughness * roughness;
        float a2 = a * a;
        float m = NdotH * NdotH * (a2 - 1) + 1;
        return a2 / std::max(M_PI * m * m, 0.000001f);
    }
    //������Ӱ�ڵ�����G
    //Disney way:
    //��G1��ʵ��
    float SmithG_GGX(float NdotV, float roughness) {
        float r = 0.5 + roughness / 2.0f;
        float m = r * r + (1 - r * r) * NdotV * NdotV;

        return 2.0f * NdotV / (NdotV + std::sqrt(m));
    }

    //��Դ����͹۲췽��ֱ����ggx1��ggx2����˵õ�G
    float GeometrySmith_Disney(Vector3f N, Vector3f V, Vector3f L, float roughness) {
        float NdotV = std::max(dotProduct(N, V), 0.0f);
        float NdotL = std::max(dotProduct(N, L), 0.0f);
        float ggx1 = SmithG_GGX(NdotL, roughness);
        float ggx2 = SmithG_GGX(NdotV, roughness);

        return ggx1 * ggx2;
    }
    //UE4 way
    //G1
    float GeometrySchlickGGX(float NdotV, float roughness) {
        float r = roughness + 1;
        float k = r * r / 8;
        float m = NdotV / NdotV * (1.f - k) + k;

        return NdotV / m;
    }

    //��Դ����͹۲췽��ֱ����ggx1��ggx2����˵õ�G
    float GeometrySmith_UE4(Vector3f N, Vector3f V, Vector3f L, float roughness) {
        float NdotV = std::max(dotProduct(N, V), 0.0f);
        float NdotL = std::max(dotProduct(N, L), 0.0f);
        float ggx1 = GeometrySchlickGGX(NdotL, roughness);
        float ggx2 = GeometrySchlickGGX(NdotV, roughness);

        return ggx1 * ggx2;
    }
    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }
        case MICROFACET:
        {
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
            return toWorld(localRay, N);
            break;
        }
        case MIRROR: 
        {
            Vector3f localRay = reflect(wi, N); 
            //�Ѿ�������������
            return localRay; 
            break;
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MICROFACET:
        {
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MIRROR: 
        {
            if (dotProduct(wo, N) > 0.0f)
                return 1.0f; 
            else
                return 0.0f;
            break; 
        }
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            //��Դ��N�ļн�
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MICROFACET:
        {
            float cosalpha = dotProduct(N, wo);//ֻ��������һ�벻��������Ҫ�ж�һ��wo��N�ļн�
            if (cosalpha > 0.0f) {
                //�����ļ���
            //ע�⣺wi������������ߵķ�����˼�����Ҫ���-wi
                float roughness = 0.8;
                Vector3f H = (-wi + wo).normalized();
                float NdotH = std::max(dotProduct(N, H), 0.0f);
                float NdotV = std::max(dotProduct(N, wo), 0.0f);
                float NdotL = std::max(dotProduct(N, -wi), 0.0f);

                //���㷨���ܶȺ��� D
                float D = DistributionGGX(NdotH, roughness);

                //������Ӱ�ڵ����� G 
                //Disney Way
                float G = GeometrySmith_Disney(N, wo, -wi, roughness);
                //UE4 way
                //float G = GeometrySmith_UE4(N, wo, -wi, roughness);

                //�������� F
                float ior = 1.5;
                float F;
                fresnel(wi, N, ior, F);

                //���㾵�淴���BRDF
                float m = 4 * std::max(NdotL * NdotL, 0.00001f);
                float Specular = D * G * F / m;

                //���㷴��������ռ��
                //����ks
                float ks = F;
                float kd = 1 - F;

                //�������BRDF
                float rou = 1.0f;//������������ʷ����ʣ���С��(0,1)��ֱ�Ӹ�1
                float Diffuse = rou / M_PI;
                //ֵ��ע����ǣ������Kd��Ks���Ƕ�Ӧ����ɫ����Specular�Ѿ��˹�F�ˣ��Ѿ������˷�����ռ�ȣ�������Ͳ����ٳ���ks��
                return Kd * Diffuse * kd + Ks * Specular;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MIRROR:
        {
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f)
            {
                float divisor = cosalpha; 
                if (divisor < 0.001) return 0;
                Vector3f mirror = 1 / divisor;
                float F;
                fresnel(wi, N, ior, F); 
                return F * mirror;
            }
            else
                return Vector3f(0.0f);
            break;
        }
    }
}

#endif //RAYTRACING_MATERIAL_H
