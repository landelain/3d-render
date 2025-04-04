#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <random>

#define pi  3.14159

double sqr(double& a ){ return a*a;};

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector operator*(const Vector& a, const Vector& b){
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}


class Ray {
public:
	Vector origin;
    Vector direction;

    Ray(const Vector& o, const Vector& d){
        origin = o;
        direction = d;
    }

};

class Sphere {
public:
    Vector origin;
    double radius;
    Vector albedo;
    bool ismirror;

	Sphere(const Vector& o, const double& r, const Vector& a, const bool& im = false){
        origin = o;
        radius = r;
        albedo = a;
        ismirror = im;
    }

    bool intersect(const Ray& ray, Vector& point, Vector& normal){
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

void get_tangent(const Vector& normal, Vector& tangent, int& idx){
    double x = normal.data[0];
    double y = normal.data[1];
    double z = normal.data[2];

    if(x < y && x < z){
        tangent.data[0] = 0.0;
        tangent.data[1] = -z;
        tangent.data[2] = y;
        idx = 0;
    }
    else if (y < x && y < z)
    {
        tangent.data[0] = -z;
        tangent.data[1] = 0.0;
        tangent.data[2] = x;
        idx = 1;
    }
    else{
        tangent.data[0] = -y;
        tangent.data[1] = x;
        tangent.data[2] = 0.0;
        idx = 2;
    }

    tangent.normalize();
}

void random_vector(const Vector& normal, Vector& random){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> unif(0.0, 1.0);
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
    int idx;
    get_tangent(normal, tangent, idx);
    tangent2 = cross(normal, tangent);
    random = coor[idx] * normal + coor[(idx + 1) % 3] * tangent + coor[(idx + 2)% 3] * tangent2;
}


class Scene {
public :
    Vector light;
    double intensity;
    std::vector<Sphere> Spheres;

    Scene(Vector& l, double& i, std::vector<Sphere>& S) {
        light = l;
        intensity = i;
        Spheres = S;
    }

    void add(Sphere& S) {
        Spheres.push_back(S);
    }

    Vector getColor(const Ray& ray, int depth, int depth2, int N) {

        if(depth2 <= 0){
            return Vector(1000, 1000, 1000);
        }

        std::vector<Sphere>::iterator it1 = Spheres.begin();
        std::vector<Sphere>::iterator it2 = Spheres.begin();

        Vector intersection;
        Vector normal;   

        bool found = false;
        Sphere BestS = *it1;
        double bestdist;

        while(it1 < Spheres.end()) {

            if(it1->intersect(ray, intersection, normal)){
                found = true;
                BestS = *it1;
                bestdist = (intersection - ray.origin).norm2();
                break;
            }
            ++it1;
        }

        if ( ! found) {
            return Vector(1000,1000,1000);
        }

        while(it2 < Spheres.end()){

            if(it2->intersect(ray, intersection, normal)){

                if((intersection - ray.origin).norm2() < bestdist){
                    BestS = *it2;
                    bestdist = (intersection - ray.origin).norm2();
                }

            }
            ++it2;
        }

        BestS.intersect(ray, intersection, normal);
        intersection = intersection + pow(10, -12)*normal;
        Vector LP = light - intersection;

        Ray shadow_ray(intersection, LP/LP.norm());
        Vector intersection2;
        Vector normal2;
        std::vector<Sphere>::iterator it3 = Spheres.begin();
        while(it3 < Spheres.end()){
            bool flag = it3->intersect(shadow_ray, intersection2, normal2);
            if(flag && (LP.norm() > (intersection2 - intersection).norm())){
                Vector color(0,0,0);
                for(int i = 0; i< N; i++){
                    Vector random;
                    random_vector(normal, random);
                    Ray rec_ray(intersection, random);
                    Vector rec_color =  getColor(rec_ray, depth, depth2-1, N);
                    rec_color = BestS.albedo * rec_color / N;
                    color = color + rec_color;
                }
                return color;
            }
            ++it3;
        }

        if(BestS.ismirror && depth > 0){
            Vector reflection_direction = ray.direction - 2*dot(ray.direction, normal) * normal;
            Ray reflection_ray(intersection, reflection_direction);
            Vector color = getColor(reflection_ray, depth-1, depth2, N);
            return color;
        }

        double dotprod = dot(normal,LP/LP.norm());
        if (dotprod < 0){return Vector(100,100,100); }
        Vector color =  intensity/(4*pi*LP.norm2()) * BestS.albedo/pi * dotprod;

        for(int i = 0; i< N; i++){
                    Vector random;
                    random_vector(normal, random);
                    Ray rec_ray(intersection, random);
                    Vector rec_color =  getColor(rec_ray, depth, depth2-1, N);
                    rec_color = BestS.albedo * rec_color / N;
                    color = color + rec_color;
                }
        return color;
    };

};


int main() {
	int W = 512;
	int H = 512;
    Vector camera(0, 0, 55);
    int max_depth = 2;
    int light_depth = 3;
    int N = 10;

    double fov = 60 * pi / 180;
    Vector albedo(1, 1, 1);
    Sphere S(Vector(0, 0, 10), 10, albedo, true);
    Sphere S2(Vector(20,10,15), 10, Vector(0,0,1));
    Vector light(-10, 20, 40);
    double intensity = 5*pow(10,9.5);

    double bigradius = 940;
    Sphere right(Vector(1000, 0, 0), bigradius, Vector(0.6, 0.5, 0.1));
    Sphere left(Vector(-1000, 0, 0), bigradius, Vector(0.9, 0.2, 0.9));
    Sphere front(Vector(0, 0, -1000), bigradius, Vector(0.4, 0.8, 0.7));
    Sphere ceil(Vector(0, 1000, 0), bigradius, Vector(0.2, 0.5, 0.9));
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.3, 0.4, 0.7));
    Sphere back(Vector(0, 0, 1000), bigradius, Vector(0.9, 0.4, 0.3), true);
    std::vector<Sphere> room;
    Scene scene(light, intensity, room);

    scene.add(ceil);
    scene.add(floor);
    scene.add(front);
    scene.add(left);
    scene.add(right);
    scene.add(back);
    scene.add(S);
    scene.add(S2);

	std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for collapse(2) schedule(guided)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            std::cout << i << std::endl;

            double z = -W/(2.*tan(fov/2));
            Vector ray_direction(j - W/2 +0.5, H/2 - i - 0.5, z);
            ray_direction.normalize();
            Ray ray(camera, ray_direction);

            Vector color = scene.getColor(ray, max_depth, light_depth, N);
            
            image[(i * W + j) * 3 + 0] = std::max(0., std::min(255., std::pow(color.data[0], 1/2.2)));
            image[(i * W + j) * 3 + 1] = std::max(0., std::min(255., std::pow(color.data[1], 1/2.2)));
            image[(i * W + j) * 3 + 2] = std::max(0., std::min(255., std::pow(color.data[2], 1/2.2)));

		}
	}
	stbi_write_png("image2.png", W, H, 3, &image[0], 0);

	return 0;
}