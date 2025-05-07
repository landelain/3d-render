#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "obj_reader.cpp"

#include <iostream>
#include <list>
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

class Properties {
public:
    bool ismirror, istrans, fresnel;

    Properties(){
        ismirror = false;
        istrans = false;
        fresnel = false;
    }

    Properties(bool im, bool it, bool fr){
        ismirror = im;
        if (im){
            istrans = false;
            fresnel = false;
            return;
        }
        istrans = it;
        fresnel = fr;
    }

};

class Object {
public:
    Vector albedo;
    bool ismirror, istrans, fresnel;


    Object(const Vector& a, const Properties& prop){
        albedo = a;
        ismirror = prop.ismirror;
        istrans = prop.istrans;
        fresnel = prop.fresnel;
    }

    virtual bool intersect(const Ray& ray, Vector& point, Vector& normal) = 0;

};

class Sphere : public Object{
public:
    Vector origin;
    double radius;

	Sphere(const Vector& o, const double& r, const Vector& a, const Properties& prop = Properties()) : Object(a, prop){
        origin = o;
        radius = r;
    }

    bool intersect(const Ray& ray, Vector& point, Vector& normal) override {
        double d = dot(ray.direction, ray.origin - origin);
        double delta = d*d - ( (ray.origin - origin).norm2() - sqr(radius));
        if(delta < 0){ return false; }
        double t1 = -d + sqrt(delta);
        double t2 = -d - sqrt(delta);
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

    void getBox(TriangleMesh& tmesh, int first, int last){
        
        max = Vector(-1e9, -1e9, -1e9);
        min = Vector(1e9, 1e9, 1e9);

        for(int i = first ; i < last; i++){
            Vector v1 = tmesh.vertices[tmesh.indices[i].vtxi];
            Vector v2 = tmesh.vertices[tmesh.indices[i].vtxj];
            Vector v3 = tmesh.vertices[tmesh.indices[i].vtxk];

            for(auto v : {v1, v2, v3}){
                max[0] = std::max(max[0], v[0]);
                max[1] = std::max(max[1], v[1]);
                max[2] = std::max(max[2], v[2]);
                
                min[0] = std::min(min[0], v[0]);
                min[1] = std::min(min[1], v[1]);
                min[2] = std::min(min[2], v[2]);
            }
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

	    if (std::min(txmax, std::min(tymax, tzmax)) > std::max(txmin, std::max(tymin, tzmin))){
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
        left = nullptr;
        right = nullptr;
        first = -1;
        last = -1;
    }

    void getBox(TriangleMesh& tmesh){
        box.getBox(tmesh, first, last);
    }
};

class NodeStack {
public:
    BVHnode* stack[100];
    int idx;

    NodeStack(){idx = -1;}

    BVHnode* pop(){
        if(idx < 0){
            return nullptr;
        }
        return stack[idx--];
    }

    void push(BVHnode* node){
        if(idx < 99){
            stack[++idx] = node;
        }
        else{
            std::cout << "stack full" << std::endl;
        }
    }

    bool isempty(){
        return idx == -1;
    }
    
};

class Mesh : public Object {
public: 
    TriangleMesh tmesh;
    BoundingBox box;
    BVHnode* root;

    Mesh(const char* l, const Vector& a, const Properties& prop = Properties()) : Object(a, prop) {
        tmesh.readOBJ(l);
        box.getBox(tmesh, 0, tmesh.indices.size());
        root = new BVHnode();
        //construct(root, 0, tmesh.indices.size());
    }

    void construct(BVHnode*& c, int first, int last){
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
                std::swap(tmesh.indices[i], tmesh.indices[pivot]);
                pivot+= 1;
            }
        }

        //std::cout<< 1.*(pivot-first)/(1*(last-first)) << std::endl;

        if (pivot <= first or pivot >= last){
            return;
        }
        if (last - first < 5){
            return;
        }

        c->left = new BVHnode();
        c->right = new BVHnode();
        construct(c->left, first, pivot);
        construct(c->right, pivot, last);
    }

    void resize(Vector& shift, double coef){
        std::vector<Vector>::iterator it = tmesh.vertices.begin();
        std::vector<Vector>::iterator end = tmesh.vertices.end();
        while(it < end){
            (*it)[0] = ( (*it)[0] - shift[0]) * coef;
            (*it)[1] = ( (*it)[1] - shift[1]) * coef;
            (*it)[2] = ((*it)[2] - shift[2]) * coef;
            it++;
        }
        box.getBox(tmesh, 0, tmesh.indices.size());
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
        NodeStack stack;
        stack.push(root);
        BVHnode* c;
        double bestt = 1e9;
        bool found = false;

        while(! stack.isempty()){
            BVHnode* c = stack.pop();
            
            if(c->left){
                if ( c->left->box.intersect(ray) ){
                    stack.push(c->left);
                }
                if ( c->right->box.intersect(ray) ){
                    stack.push(c->right);
                }
            } else {
                for (int i = c->first; i < c->last; i++) {
                    const Vector& A = tmesh.vertices[tmesh.indices[i].vtxi];
                    const Vector& B = tmesh.vertices[tmesh.indices[i].vtxj];
                    const Vector& C = tmesh.vertices[tmesh.indices[i].vtxk];
                    double alpha, beta, gamma, t;

                    if (MollerTrumbore(ray, A, B, C, t, alpha, beta, gamma)){
                        if(t < bestt){
                            found = true;
                            bestt = t;
                            point = alpha * A + beta * B + gamma * C;
                            normal = alpha*tmesh.normals[tmesh.indices[i].ni] + beta*tmesh.normals[tmesh.indices[i].nj] + gamma*tmesh.normals[tmesh.indices[i].nk];
                            normal.normalize();
                        }
                        
                    } 
                }
            }  
        }

        return found;
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

    Vector getColor(const Ray& ray, int depth) {

        Vector intersection;
        Vector normal;   

        bool found = false;
        Object* BestO = Objects[0];
        double bestdist = 1e9;
        Vector bestinter, bestnorm;

        for( auto O : Objects){
            if((O->intersect(ray, intersection, normal))){
                    found = true;
                if ((intersection - ray.origin).norm2() < bestdist){
                    BestO = O;
                    bestdist = (intersection - ray.origin).norm2();
                    bestinter = intersection;
                    bestnorm = normal;
                }
            }
        }

        if ( ! found) {
            return Vector(100,100,100);
        }

        double eps = pow(10, -3);
        intersection = bestinter;
        normal = bestnorm;
        Vector LP = light - intersection;

        if(BestO->ismirror && depth > 0){
            Vector reflection_direction = ray.direction - 2*dot(ray.direction, normal) * normal;
            Ray reflection_ray(intersection + eps * reflection_direction, reflection_direction);
            Vector color = getColor(reflection_ray, depth-1);
            return color;
        }
        else if(BestO->istrans && depth > 0){
            double n1 = 1;
            double n2 = 1.5;

            // Snell
            double iN = dot(ray.direction, normal);
            if(iN > 0){
                // Inside the sphre so swap the mediums / normal
                std::swap(n1, n2);
                normal = -1. * normal;
                iN = dot(ray.direction, normal);
            }
            Vector trans_direction;

            double p = 1-(n1/n2)*(n1/n2)*(1-iN*iN);
            if(p < 0){
                // Since negative we get full reflection like a mirror
                trans_direction = ray.direction - 2*iN * normal;
                trans_direction.normalize();
                Ray trans_ray(intersection + eps * trans_direction, trans_direction);
                Vector color = getColor(trans_ray, depth-1);
                return color;
            }
            Vector totT = n1/n2 * (ray.direction - iN * normal);
            Vector tonN =  - sqrt(p) * normal;
            trans_direction = totT + tonN;

            // Fresnel
            if(BestO->fresnel){
                double k0 = pow(n1-n2, 2)/pow(n1 + n2, 2);
                double R = k0 + (1 - k0)*pow(1 - abs(iN), 5);
                if(unif(gen) < R){
                    trans_direction = ray.direction - 2*iN * normal;
                }
            }

            trans_direction.normalize();
            Ray trans_ray(intersection + eps * trans_direction, trans_direction);
            Vector color = getColor(trans_ray, depth-1);
            return color;
        }

        intersection = intersection + eps*normal;
        Ray shadow_ray(intersection, LP/LP.norm());
        Vector intersection2;
        Vector normal2;
        bool shadow = false;
        double dotprod = dot(normal, LP/LP.norm()); 
        if (dotprod < 0) {
            shadow = true;
        }
        else{
            for(auto O : Objects){
                bool flag = O->intersect(shadow_ray, intersection2, normal2);
                if(flag && (LP.norm2() > (intersection2 - intersection).norm2())){
                    shadow = true;
                    break;
                }
            }
        }

        Vector color(0,0,0);
        if (! shadow){ 
            double dotprod = dot(normal,LP/LP.norm());
            if (dotprod < 0){return Vector(0,0,0); }
            color =  intensity/(4*pi*LP.norm2()) * BestO->albedo/pi * dotprod;
        }

        if(depth <= 0 ){return color;}

        Vector random;
        random_vector(normal, random);
        Ray rec_ray(intersection, random);
        color =  color + BestO->albedo * getColor(rec_ray, depth-1);
    
        return color;
    };

};


int main() {
	int W = 512;
	int H = 512;
    // int W = 256;
    // int H = 256;

    Vector camera(0, 0, 55);
    int depth = 5;
    int N = 64;

    double fov = 60 * pi / 180;
    Vector albedo(1, 0, 0);
    Properties prop(true, false, false);
    Properties prop2(false, true, true);
    Sphere S(Vector(0, 0, 10), 7, albedo, prop);
    Sphere S2(Vector(-15,0,10), 7, albedo);
    Sphere S3(Vector(15,0,10), 7, albedo, prop2);
    Vector light(-10, 20, 40);
    double intensity = 5*pow(10,9);

    double bigradius = 940;
    Sphere right(Vector(1000, 0, 0), bigradius, Vector(0.6, 0.5, 0.1));
    Sphere left(Vector(-1000, 0, 0), bigradius, Vector(0.9, 0.2, 0.9));
    Sphere front(Vector(0, 0, -1000), bigradius, Vector(0.4, 0.8, 0.7));
    Sphere ceil(Vector(0, 1000, 0), bigradius, Vector(0.2, 0.5, 0.9));
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.3, 0.4, 0.7));
    Sphere back(Vector(0, 0, 1000), bigradius, Vector(0.9, 0.4, 0.3));
    std::vector<Object*> room;
    Scene scene(light, intensity, room);

    Properties propcat(false, false, false);
    Mesh cat_mesh = Mesh("cat_model/cat.obj", Vector(0.8,0.8,0.8), propcat);
    Vector resize(0, 15, 0);
    cat_mesh.resize(resize, 0.4);
    cat_mesh.construct(cat_mesh.root, 0, cat_mesh.tmesh.indices.size());
    scene.add(&cat_mesh);
    
    scene.add(&ceil);
    scene.add(&floor);
    scene.add(&front);
    scene.add(&left);
    scene.add(&right);
    scene.add(&back);
    //scene.add(&S);
    // scene.add(&S2);
    // scene.add(&S3);

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
                color = color + scene.getColor(ray, depth);
            }
            color = color / N;
            
            image[(i * W + j) * 3 + 0] = std::max(0., std::min(255., std::pow(color.data[0], 1./2.2)));
            image[(i * W + j) * 3 + 1] = std::max(0., std::min(255., std::pow(color.data[1], 1./2.2)));
            image[(i * W + j) * 3 + 2] = std::max(0., std::min(255., std::pow(color.data[2], 1./2.2)));

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

    //delete cat_mesh;
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
