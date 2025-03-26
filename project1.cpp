#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

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

	Sphere(const Vector& o, const double& r, const Vector& a){
        origin = o;
        radius = r;
        albedo = a;
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

class Scene {
public :
    std::vector<Sphere> Spheres;

    Scene(std::vector<Sphere>& S) {
        Spheres = S;
    }

    void add(Sphere& S) {
        Spheres.push_back(S);
    }

};


Vector getColor(Vector& light, double& intensity, Vector& albedo, Vector& normal, Vector& point) {
    Vector LP = light - point;
    double dotprod = dot(normal,LP/LP.norm());
    if (dotprod < 0){return Vector(0,0,0); }
    return intensity/(4*pi*LP.norm2()) * albedo/pi * dotprod;
};


int main() {
	int W = 512;
	int H = 512;
    Vector camera(0, 0, 55);

    double fov = 60 * pi / 180;
    Vector albedo(1, 0, 0);
    Sphere S(Vector(0, 0, 0), 10, albedo);
    Vector light(20, 20, 20);
    double intensity = 10000000;

    Sphere ceil(Vector(1000, 0, 0), 940, Vector(0.2, 0.5, 0.9));
    Sphere floor(Vector(-1000, 0, 0), 940, Vector(0.3, 0.4, 0.7));
    Sphere front(Vector(0, 0, -1000), 940, Vector(0.4, 0.8, 0.7));
    Sphere left(Vector(0, 1000, 0), 940, Vector(0.9, 0.2, 0.9));
    Sphere right(Vector(0, -1000, 0), 940, Vector(0.6, 0.5, 0.1));
    Sphere back(Vector(0, 0, 1000), 940, Vector(0.9, 0.4, 0.3));
    std::vector<Sphere> room;
    Scene scene(room);

    //scene.add(ceil);
    //scene.add(floor);
    //scene.add(front);
    //scene.add(left);
    //scene.add(right);
    //scene.add(back);
    scene.add(S);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            double z = -W/(2.*tan(fov/2));
            Vector ray_direction(j - W/2 +0.5, H/2 - i - 0.5, z);
            ray_direction.normalize();
            Ray ray(camera, ray_direction);
            Vector intersection(0,0,0);
            Vector normal(0,0,0);

            bool colored = false;
            std::vector<Sphere>::iterator it = scene.Spheres.begin();

            while(it < scene.Spheres.end()) {

                if(it->intersect(ray, intersection, normal)){
                    colored = true;
                    
                    Vector color = getColor(light, intensity, albedo, normal, intersection);           

                    image[(i * W + j) * 3 + 0] = std::max(0., std::min(255., color.data[0]));
                    image[(i * W + j) * 3 + 1] = std::max(0., std::min(255., color.data[1]));
                    image[(i * W + j) * 3 + 2] = std::max(0., std::min(255., color.data[2]));

                    break;
                }

                ++it;
            }

            if(! colored) {
                image[(i * W + j) * 3 + 0] = 100;
                image[(i * W + j) * 3 + 1] = 100;
                image[(i * W + j) * 3 + 2] = 100;
            }
			
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}