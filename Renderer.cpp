//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include <atomic>
#include <mutex>
#include <thread>
#include <omp.h>

#pragma warning(disable : 4996)

inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

std::mutex mtx;
const float EPSILON = 0.001;
// ����̶߳�ͬһ�������޸���Ҫ�������� ������������̫�˷�ʱ�� ����atomicͬ������
//std::atomic_int progress = 0;


// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;
    
    // change the spp value to change sample ammount
    int spp = 512;
    std::cout << "SPP: " << spp << "\n";

    int progress = 0;
    int thred = 24;
    // һ���̴߳���40��
    int times = scene.height / thred;  // 960/24=40
    std::thread th[24]; //���߳�

    int width, height;
    width = height = sqrt(spp);
    float step = 1.0f / width;

    // �߳�ʹ���� C++11��׼������lambda����
    // [&]��ʾ����ǰ���������б��������ô���
    
    auto castRayMultiThread = [&](uint32_t lrow, uint32_t hrow) {
    //#pragma omp parallel for
        for (int j = lrow; j < hrow; ++j) {
            //int m = j * scene.width;
            for (uint32_t i = 0; i < scene.width; ++i) {
                // generate primary ray direction
                /*float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                    imageAspectRatio * scale;
                float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

                Vector3f dir = normalize(Vector3f(-x, y, 1));*/
                for (int k = 0; k < spp; k++) {
                    float x = (2 * (i + step /2+ step *(k% width)) / (float)scene.width - 1) * imageAspectRatio * scale;
                    float y = (1 - 2 * (j + step /2 + step *(k/height)) / (float)scene.height) * scale;

                    Vector3f dir = normalize(Vector3f(-x, y, 1));
                    framebuffer[(int)(j * scene.width + i)] += scene.castRay(Ray(eye_pos, dir),0) / spp;
                }
            }
            mtx.lock();
            progress++;
            UpdateProgress(progress / (float)scene.height);
            mtx.unlock();
        }
    };
    //���н���·��׷��
    //i*times��ʾ���ǵ�i���߳�������ʼֵ��(i+1)*times��ʾ��i���߳�����������У��ֱ��Ӧ������y_min��y_max.
    #pragma omp parallel for
    for (int i = 0; i < thred; i++) {//�ӵ�0�г�����һ����0~by-1��
        th[i] = std::thread(castRayMultiThread, i * times, (i + 1) * times);
    }
    //û��ִ��join��detach���߳��ڳ������ʱ�������쳣�������Ҫ��ÿ���̶߳�ִ��һ��join.
    //ÿ���߳�ִ��join
    #pragma omp parallel for
    for (int i = 0; i < thred; i++) {
        th[i].join();
    }
    
    
    // change the spp value to change sample ammount
    //int spp = 16;
    //std::cout << "SPP: " << spp << "\n";
    //#pragma omp parallel for
    //for (int j = 0; j < scene.height; ++j) {
    //    for (uint32_t i = 0; i < scene.width; ++i) {
    //        // generate primary ray direction
    //        float x = (2 * (i + 0.5) / (float)scene.width - 1) *
    //            imageAspectRatio * scale;
    //        float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

    //        Vector3f dir = normalize(Vector3f(-x, y, 1));
    //        for (int k = 0; k < spp; k++) {
    //            framebuffer[(int)(j * scene.width + i)] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
    //            //framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
    //        }
    //        //m++;
    //    }
    //    UpdateProgress(j / (float)scene.height);
    //}
    
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}