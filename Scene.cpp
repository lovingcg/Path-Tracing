//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    // 发光区域面积
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum) {
                //随机选取一个光源面，即第k个自发光物体的光源面
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

/*
// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Vector3f L_dir;
    Vector3f L_indir;

    // 从像素发出的光线与物体的交点
    Intersection obj_inter = intersect(ray);
    if (!obj_inter.happened)
        return L_dir;

    // 打到光源
    if (obj_inter.m->hasEmission())
        return obj_inter.m->getEmission();

    // 打到物体
    Vector3f p = obj_inter.coords;
    Material* m = obj_inter.m;
    Vector3f N = obj_inter.normal.normalized();
    Vector3f wo = ray.direction; // 像素到物体的向量

    // 有交点，对光源采样

    float pdf_L = 1.0; //可以不初始化
    Intersection light_inter;
    sampleLight(light_inter, pdf_L);    // 得到光源位置和对光源采样的pdf

    Vector3f x = light_inter.coords;
    Vector3f ws = (x - p).normalized(); //物体到光源

    Vector3f NN = light_inter.normal.normalized();
    Vector3f emit = light_inter.emit;
    float d = (x - p).norm();

    // 再次从光源发出一条光线，判断是否能打到该物体，即中间是否有阻挡

    Ray Obj2Light(p, ws);
    float d2 = intersect(Obj2Light).distance;
    // 是否阻挡，利用距离判断，需注意浮点数的处理
    // 未阻挡 直接光照
    if (d2 - d > -0.001) {
        Vector3f eval = m->eval(wo, ws, N); // wo不会用到
        float cos_theta = dotProduct(N, ws);
        float cos_theta_x = dotProduct(NN, -ws);//ws从物体指向光源，与NN的夹角大于180
        // emit表示Li(p,wi) eval表示BRDF cos_theta表示n点乘wi cos_theta_x表示n撇点乘wi
        L_dir = emit * eval * cos_theta * cos_theta_x / std::pow(d, 2) / pdf_L;
    }

    // L_indir
    float P_RR = get_random_float();
    if (P_RR < RussianRoulette) {
        //随机生成一个wi方向
        Vector3f wi = m->sample(wo, N).normalized();
        Ray r(p, wi);
        Intersection inter = intersect(r);
        // 判断打到的物体是否会发光取决于m
        if (inter.happened && !inter.m->hasEmission()) {
            Vector3f eval = m->eval(wo, wi, N);
            float pdf_O = m->pdf(wo, wi, N);
            float cos_theta = dotProduct(wi, N);
            // 对应间接光照的渲染公式
            L_indir = castRay(r, depth + 1) * eval * cos_theta / pdf_O / RussianRoulette;
        }
    }
    //4->16min
    return L_dir + L_indir;
}
*/

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray& ray, int depth) const
{
    //创建变量以储存直接和间接光照计算值

    Vector3f dir = { 0.0,0.0,0.0 };
    Vector3f indir = { 0.0,0.0,0.0 };

    //1.判断是否有交点：光线与场景中物体相交？

    Intersection inter = Scene::intersect(ray);
    //如果没交点

    if (!inter.happened) {
        return dir;//return 0,0,0
    }
    //2.ray打到光源了：说明渲染方程只用算前面的自发光项，因此直接返回材质的自发光项 

    // 打到光源
    if (inter.m->hasEmission())
        return inter.m->getEmission();

    //3.ray打到物体：这个时候才开始进行伪代码后面的步骤

    //对场景中的光源进行采样，得到采样点light_pos和pdf_light
    Intersection light_pos;
    float pdf_light = 0.0f;
    sampleLight(light_pos, pdf_light);

    //3.1计算直接光照
    //物体的一些参数

    Vector3f p = inter.coords;
    Vector3f N = inter.normal.normalized();
    Vector3f wo = ray.direction;//物体指向场景
    
    switch (inter.m->getType())
    {
        case MIRROR:
        {
            float ksi = get_random_float();//随机取[0,1]
            if (ksi < RussianRoulette) {
                //随机生成一个wi方向
                Vector3f wi = inter.m->sample(wo, N).normalized();//这里的wi其实没参与计算，返回的是一个随机的方向
                Ray r(p, wi);
                Intersection obj_to_scene = Scene::intersect(r);
                if (obj_to_scene.happened)
                {
                    Vector3f f_r = inter.m->eval(wo, wi, N);//wo不参与计算

                    float cos_theta = dotProduct(wi, N);
                    //float pdf_hemi = inter.m->pdf(wo, wi, N);
                    float pdf_hemi = inter.m->pdf(wo, wi, N);
                    if (pdf_hemi > EPSILON)
                    {
                        indir = castRay(r, depth + 1) * f_r * cos_theta / pdf_hemi / RussianRoulette;
                    }

                }
            }
            break;
        }
        default:
        {
            Vector3f xx = light_pos.coords;
            Vector3f NN = light_pos.normal.normalized();
            Vector3f ws = (p - xx).normalized();//光源指向物体
            float dis = (p - xx).norm();//二者距离

            float dis2 = dotProduct((p - xx), (p - xx));

            //判断光源与物体间是否有遮挡：
            //发出一条射线，方向为ws 光源xx -> 物体p

            Ray light_to_obj(xx, ws);//Ray(orig,dir)
            Intersection light_to_scene = Scene::intersect(light_to_obj);
            //假如dis>light_to_scene.distance就说明有遮挡，那么反着给条件即可：
            if (light_to_scene.happened && (light_to_scene.distance - dis >= -sqrt(EPSILON))) {//没有遮挡
                //为了更贴近伪代码，先设定一些参数

                Vector3f L_i = light_pos.emit;//光强
                Vector3f f_r = inter.m->eval(wo, -ws, N);//材质，课上说了，BRDF==材质，ws不参与计算

                float cos_theta = dotProduct(-ws, N);//物体夹角
                float cos_theta_l = dotProduct(ws, NN);//光源夹角
                dir = L_i * f_r * cos_theta * cos_theta_l / dis2 / pdf_light;
            }

            float ksi = get_random_float();//随机取[0,1]
            if (ksi < RussianRoulette) {
                //计算间接光照

                //随机生成一个wi方向
                Vector3f wi = inter.m->sample(wo, N).normalized();//这里的wi其实没参与计算，返回的是一个随机的方向
                Ray r(p, wi);
                Intersection obj_to_scene = Scene::intersect(r);
                if (obj_to_scene.happened && !obj_to_scene.m->hasEmission())
                {
                    Vector3f f_r = inter.m->eval(wo, wi, N);//wo不参与计算

                    float cos_theta = dotProduct(wi, N);
                    float pdf_hemi = inter.m->pdf(wo, wi, N);
                    if (pdf_hemi > EPSILON)
                    {
                        indir = castRay(r, depth + 1) * f_r * cos_theta / pdf_hemi / RussianRoulette;
                    }
                }
            }
            break;
        }
    }
    Vector3f final_color= dir + indir;
    final_color.x = (clamp(0, 1, final_color.x));
    final_color.y = (clamp(0, 1, final_color.y));
    final_color.z = (clamp(0, 1, final_color.z));
    return final_color;
    //return dir + indir;
}
