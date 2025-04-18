#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "obj_reader.cpp"

#include <iostream>
#include <random>

#define pi  3.14159

double sqr(double a ){ return a*a;};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unif(0.0, 1.0);

class Ray {
public:
	Vector origin;
    Vector direction;

    Ray(const Vector& o, const Vector& d){
        origin = o;
        direction = d;
    }

};

class Object {
public:
    Vector albedo;
    bool ismirror;

    Object(const Vector& a, const bool& im = false){
        albedo = a;
        ismirror = im;
    }

    virtual bool intersect(const Ray& ray, Vector& point, Vector& normal) = 0;

};

class Sphere : public Object{
public:
    Vector origin;
    double radius;

	Sphere(const Vector& o, const double& r, const Vector& a, const bool& im = false) : Object(a, im){
        origin = o;
        radius = r;
    }

    bool intersect(const Ray& ray, Vector& point, Vector& normal) override {
        double delta = (dot(ray.direction, ray.origin - origin))*(dot(ray.direction, ray.origin - origin)) - ( (ray.origin - origin).norm2() - sqr(radius));
        if(delta < 0){ return false; }
        double dotprod = dot(ray.direction, origin - ray.origin);
        double t1 = dotprod + sqrt(delta);
        double t2 = dotprod - sqrt(delta);
        if (t1 < 0){ return false;}
        double t;
        if( t2 > 0 ){
            t = t2;
        }
        else{
            t = t1;
        }
        point = ray.origin + t*ray.direction;
        normal = point - origin;
        normal.normalize();
        return true;
    }
};

class BoundingBox {
public:

    Vector max, min;

    BoundingBox(){}

    BoundingBox(Vector& max, Vector& min){
        this->max = max;
        this->min = min;
    }

    void getBox(TriangleMesh& tmesh, int first, int last){
        max = Vector(-1e9, -1e9, -1e9);
        min = Vector(1e9, 1e9, 1e9);

        std::vector<TriangleIndices>::iterator it1 = tmesh.indices.begin();
        std::vector<TriangleIndices>::iterator end = tmesh.indices.end();
        while( it1 < end){
            max[0] = std::max(max[0], (*it1)[0]);
            max[1] = std::max(max[1], (*it1)[1]);
            max[2] = std::max(max[2], (*it1)[2]);
            
            min[0] = std::min(min[0], (*it1)[0]);
            min[1] = std::min(min[1], (*it1)[1]);
            min[2] = std::min(min[2], (*it1)[2]);

            ++it1;
        }
    }

    bool intersect(const Ray& ray){
        double tx1 = (min[0] - ray.origin[0])/ray.direction[0];
        double tx2 = (max[0] - ray.origin[0])/ray.direction[0];
        double txmin = std::min(tx1, tx2);
        double txmax = std::max(tx1, tx2);

        double ty1 = (min[1] - ray.origin[1])/ray.direction[1];
        double ty2 = (max[1] - ray.origin[1])/ray.direction[1];
        double tymin = std::min(ty1, ty2);
        double tymax = std::max(ty1, ty2);

        double tz1 = (min[2] - ray.origin[2])/ray.direction[2];
        double tz2 = (max[2] - ray.origin[2])/ray.direction[2];
        double tzmin = std::min(tz1, tz2);
        double tzmax = std::max(tz1, tz2);

        double tmax = std::min(txmax, tymax);
        double minmax = std::min(tmax, tzmax);
        double tmin = std::max(txmin, tymin);
        double maxmin = std::max(tmin, tzmin);
        if(minmax > maxmin){
            return true;
        }
        return false;
    }
};

class BVHnode {
public:
    BoundingBox box;
    BVHnode *left, *right;
    int first, last;

    BVHnode(){

    }
    BVHnode(BVHnode* c, int first, int last){
        this->first = first;
        this->last = last;
    }

    void getBox(TriangleMesh& tmesh){
        box.getBox(tmesh, first, last);
    }
};


class Mesh : public Object {
public: 
    TriangleMesh tmesh;
    BoundingBox box;
    BVHnode* root;

    Mesh(const char* l, const Vector& a, const bool& im = false) : Object(a, im) {
        tmesh.readOBJ(l);
        box.getBox(tmesh, 0, tmesh.indices.size());
        construct(root, 0, tmesh.indices.size());
        std::cout << "passed" << std::endl;
    }

    void resize(Vector& shift, double coef){
        std::vector<Vector>::iterator it = tmesh.vertices.begin();
        while(it < tmesh.vertices.end()){

            (*it)[0] = ( (*it)[0] - shift[0]) * coef;
            (*it)[1] = ( (*it)[1] - shift[1]) * coef;
            (*it)[2] = ( (*it)[2] - shift[2]) * coef;
            ++it;
        }
        box.getBox(tmesh, 0, tmesh.indices.size());
    }


    void construct(BVHnode* c, int first, int last){
    
    //std::cout << last-first << std::endl;

    if(last-first <= (tmesh.indices.size()*0.3)){
        c = nullptr;
        return;
    }

    c = new BVHnode();

    c->first = first;
    c->last = last;
    c->getBox(tmesh);
    
    int axis = 0;
    for(int i = 1; i<3; i++){
        if(std::abs(c->box.max[i] - c->box.min[i]) > std::abs(c->box.max[axis] - c->box.min[axis])){
            axis = i;
        }
    }

    double M = (c->box.max[axis] + c->box.min[axis])/2;
    int pivot = first;
    for(int i = first; i < last; i++){
        Vector barycenter = (tmesh.vertices[tmesh.indices[i].vtxi] + tmesh.vertices[tmesh.indices[i].vtxj] + tmesh.vertices[tmesh.indices[i].vtxk])/3;
        if( barycenter[axis] < M){
            std::cout << "happens" << std::endl;
            TriangleIndices temp = tmesh.indices[i];
            tmesh.indices[i] = tmesh.indices[pivot];
            tmesh.indices[pivot] = temp;
            pivot++;
        }
    }
    
    construct(c->left, first, pivot);
    construct(c->right, pivot, last);
}

    bool MollerTrumbore(const Ray& ray, const Vector& A, const Vector& B, const Vector& C, double& t, double& alpha, double& beta, double& gamma){
        Vector O = ray.origin;
        Vector u = ray.direction;
        Vector e1 = B - A;
        Vector e2 = C - A;
        Vector N = cross(e1, e2);
        Vector AOu = cross((A - O), u);
        double eps = 1e-8;
        double uN = dot(u, N);
        //if(uN < eps){return false;}

        beta = dot(e2, AOu)/uN;
        if(beta > 1 + eps || beta < -eps){
            return false;
        }
        gamma = -dot(e1, AOu)/uN;
        if(gamma > 1 + eps || gamma < -eps){
            return false;
        }
        alpha = 1 - beta - gamma;
        if(alpha > 1 + eps || alpha < -eps){
            return false;
        }
        t = dot((A-O), N)/uN;
        if(t < 0){
            return false;
        }
        return true;
    }

    bool intersect(const Ray& ray, Vector& point, Vector& normal) override {

        // intersection with the big Bounding Box
        if(! box.intersect(ray)){
            return false;
        }

        // intersection with the node Boxes
        printf("hey");
        BVHnode* nodes[100];
        int idx = 0;
        nodes[idx] = root;

        while(idx >= 0){
            BVHnode* c = nodes[idx];
            idx--;
            if(c->left != nullptr){
                if(c->left->box.intersect(ray)){
                    idx++;
                    nodes[idx] = c->left;
                }
            }
            if(c->right != nullptr)
            {
                if(c->right->box.intersect(ray)){
                    idx++;
                    nodes[idx] = c->right;
                }
            }
            if(c->left == nullptr && c->right == nullptr)
            {
                break;
            }
            std::cout << idx << std::endl;
        }   

        int first = nodes[idx+1]->first;
        int last = nodes[idx+1]->last;

        double bestt = 10000;
        double alpha, beta, gamma, t;
        Vector bestNormal;
        bool found = false;
        Vector bestPoint;

        for(int i = first; i < last; i++){
            TriangleIndices ti = tmesh.indices[i];
            Vector A = tmesh.vertices[ti.vtxi];
            Vector B = tmesh.vertices[ti.vtxj];
            Vector C = tmesh.vertices[ti.vtxk];

            if(MollerTrumbore(ray, A, B, C, t, alpha, beta, gamma)){
                found = true;
                if (t < bestt){
                    bestt = t;
                    bestPoint = alpha * A + beta * B + gamma * C;
                    bestNormal = alpha * tmesh.normals[ti.ni] + beta * tmesh.normals[ti.nj] + gamma * tmesh.normals[ti.nk];
                } 
            }
        }  

        if(! found){
            return false;
        }

        point = bestPoint;
        normal = bestNormal;
        normal.normalize();     
        return true;
    }

};

void get_tangent(const Vector& normal, Vector& tangent){
    double x = normal.data[0];
    double y = normal.data[1];
    double z = normal.data[2];

    double ax = abs(x);
    double ay = abs(y);
    double az = abs(z);

    if(ax < ay && ax < az){
        tangent.data[0] = 0.0;
        tangent.data[1] = -z;
        tangent.data[2] = y;
    }
    else if (ay < ax && ay < az)
    {
        tangent.data[0] = -z;
        tangent.data[1] = 0.0;
        tangent.data[2] = x;
    }
    else{
        tangent.data[0] = -y;
        tangent.data[1] = x;
        tangent.data[2] = 0.0;
    }

    tangent.normalize();
}

void random_vector(const Vector& normal, Vector& random){
    double r1 = unif(gen);
    double r2 = unif(gen);

    double x = cos(2*pi*r1)*sqrt(1-r2);
    double y = sin(2*pi*r1)*sqrt(1-r2);
    double z = sqrt(r2);
    std::vector<float> coor(3);
    coor[0] = x;
    coor[1] = y;
    coor[2] = z;

    Vector tangent, tangent2;
    get_tangent(normal, tangent);
    tangent2 = cross(normal, tangent);
    random = coor[2] * normal + coor[0] * tangent + coor[1] * tangent2;
}

class Scene {
public :
    Vector light;
    double intensity;
    std::vector<Object*> Objects;

    Scene(Vector& l, double& i, std::vector<Object*>& Ob) {
        light = l;
        intensity = i;
        Objects = Ob;
    }

    void add(Object* Ob) {
        Objects.push_back(Ob);
    }

    Vector getColor(const Ray& ray, int depth, int depth2) {

        Vector intersection;
        Vector normal;   

        bool found = false;
        Object* BestO = Objects[0];
        double bestdist = 100000;

        std::vector<Object*>::iterator it1 = Objects.begin();
        while( it1 < Objects.end()){
            if(((*it1)->intersect(ray, intersection, normal))){
                    found = true;
                if ((intersection - ray.origin).norm() < bestdist){
                    BestO = *it1;
                    bestdist = (intersection - ray.origin).norm();
                }
            }
            ++it1;
        }

        if ( ! found) {
            return Vector(100,100,100);
        }

        BestO->intersect(ray, intersection, normal);
        intersection = intersection + pow(10, -4)*normal;
        Vector LP = light - intersection;

        if(BestO->ismirror && depth > 0){
            Vector reflection_direction = ray.direction - 2*dot(ray.direction, normal) * normal;
            Ray reflection_ray(intersection, reflection_direction);
            Vector color = getColor(reflection_ray, depth-1, depth2);
            return color;
        }

        Ray shadow_ray(intersection, LP/LP.norm());
        Vector intersection2;
        Vector normal2;
        bool shadow = false;
        std::vector<Object*>::iterator it2 = Objects.begin();
        while(it2 < Objects.end()){
            bool flag = (*it2)->intersect(shadow_ray, intersection2, normal2);
            if(flag && (LP.norm() > (intersection2 - intersection).norm())){
                shadow = true;
                break;
            }
            ++it2;
        }

        Vector color(0,0,0);
        if (! shadow){
            double dotprod = dot(normal,LP/LP.norm());
            if (dotprod < 0){return Vector(0,0,0); }
            color =  intensity/(4*pi*LP.norm2()) * BestO->albedo/pi * dotprod;
        }

        if(depth2 <= 0 ){return color;}

        Vector random;
        random_vector(normal, random);
        Ray rec_ray(intersection, random);
        color =  color + BestO->albedo * getColor(rec_ray, depth, depth2-1);
    
        return color;
    };

};


int main() {
	//int W = 512;
	//int H = 512;
    int W = 256;
    int H = 256;

    Vector camera(0, 0, 55);
    int max_depth = 2;
    int light_depth = 3;
    int N = 10;

    double fov = 60 * pi / 180;
    Vector albedo(1, 0, 0);
    Sphere S(Vector(0, 0, 10), 10, albedo, true);
    Sphere S2(Vector(20,10,15), 10, Vector(0,0,1));
    Vector light(-10, 20, 40);
    double intensity = 5*pow(10,9);

    double bigradius = 940;
    Sphere right(Vector(1000, 0, 0), bigradius, Vector(0.6, 0.5, 0.1));
    Sphere left(Vector(-1000, 0, 0), bigradius, Vector(0.9, 0.2, 0.9));
    Sphere front(Vector(0, 0, -1000), bigradius, Vector(0.4, 0.8, 0.7));
    Sphere ceil(Vector(0, 1000, 0), bigradius, Vector(0.2, 0.5, 0.9));
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.3, 0.4, 0.7));
    Sphere back(Vector(0, 0, 1000), bigradius, Vector(0.9, 0.4, 0.3), true);
    std::vector<Object*> room;
    Scene scene(light, intensity, room);

    Mesh* cat_mesh = new Mesh("cat_model/cat.obj", Vector(0.8,0.8,0.8));
    Vector resize(0, 15, 0);
    cat_mesh->resize(resize, 0.4);
    scene.add(cat_mesh);
    
    scene.add(&ceil);
    scene.add(&floor);
    scene.add(&front);
    scene.add(&left);
    scene.add(&right);
    scene.add(&back);
    //scene.add(&S);
    //scene.add(S2);

	std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for collapse(2) schedule(guided)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vector color(0,0,0);
            double z = -W/(2.*tan(fov/2));
            for(int n = 0; n < N; ++n){
                double r1 = unif(gen);
                double r2 = unif(gen);
                double x2 = sqrt(-2*log(r1)) * cos(2*pi*r1) * 0.25;
                double y2 = sqrt(-2*log(r1)) * sin(2*pi*r1) * 0.25;

                Vector ray_direction(j - W/2 +0.5 + x2, H/2 - i - 0.5 + y2, z);
                ray_direction.normalize();
                Ray ray(camera, ray_direction);
                color = color + scene.getColor(ray, max_depth, light_depth);
            }
            color = color / N;
            
            image[(i * W + j) * 3 + 0] = std::max(0., std::min(255., std::pow(color.data[0], 1./2.2)));
            image[(i * W + j) * 3 + 1] = std::max(0., std::min(255., std::pow(color.data[1], 1./2.2)));
            image[(i * W + j) * 3 + 2] = std::max(0., std::min(255., std::pow(color.data[2], 1./2.2)));

		}
	}
	stbi_write_png("image4.png", W, H, 3, &image[0], 0);

    delete cat_mesh;
	return 0;
}

// std vector std vector double textures
// std vector int texture W
//  textureH

// //load texture (char* file)
// int W? H, C
// unsigner char * texture = stbi_load(file, &W, &H, &C, 3)
// textureW.pushback(W);
// textureH.pushback
// std vector double current_tex(W*H*3)
// for i < W*H*3
// curretn_text[i] = std::power(tex[i] /  255. , 2.2)
// textures.push_back(current_text)

// mesh.load_texture(filename)

// return current albdedo at intersetc and not just albedo of the cat itself
// vector UV = interpolation of the UV coordiantes (much like normal)
// U = std::min(UV[0]*texturesW[indices[i].group] -1. , UV[0]*texturesW[indices[i].group]); same min for the height after and use max than 0
// V = (1-UV[1])*texturesH[indices[i].group]; maybe 1- maybe not
// texH = texturesH[indices[i].group]
// texW = same
// albedo = Vector( textures[indicies[i].group][V*texH + U + 0], textures[indicies[i].group][V*texH + U + 1], textures[indicies[i].group][V*texH + U + 2]

// dont forget that Mesh returns this albedo so OBject should also and so should Sphere




