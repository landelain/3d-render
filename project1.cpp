#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

	Sphere(const Vector& o, const double& r){
        origin = o;
        radius = r;
    }

    bool intersect(const Ray& ray, Vector& point, Vector& normal){
        double delta = sqrt(dot(ray.direction, ray.origin -origin)) - ( (ray.origin - origin).norm2() - sqr(radius));
        if(delta < 0){ return false; }
        double dotprod = dot(ray.direction, origin - ray.origin);
        double t1 = dotprod + sqrt(delta);
        double t2 = dotprod - sqrt(delta);
        if (t2 < 0){ return false;}
        double t;
        if( t1 > 0 ){
            t = t1;
        }
        else{
            t = t2;
        }
        point = ray.origin + t*ray.direction;
        normal = point - origin;
        return true;
    }
};

int main() {
	int W = 512;
	int H = 512;
    Vector camera(0, 0, 50);

    double fov = 60 * pi / 180;
    Sphere S(Vector(0, 0, 0), 10);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            double d = -W/(2.*tan(fov/2));
            Vector ray_direction(j - W/2 +0.5, H/2 - i - 0.5, d);
            ray_direction.normalize();
            Ray ray(camera, ray_direction);
            Vector intersection(0,0,0);
            Vector normal(0,0,0);

            if(S.intersect(ray, intersection, normal)){
                image[(i * W + j) * 3 + 0] = 0;
			    image[(i * W + j) * 3 + 1] = 0;
			    image[(i * W + j) * 3 + 2] = 0;
            }
            else{
                image[(i * W + j) * 3 + 0] = 255;
                image[(i * W + j) * 3 + 1] = 255;
                image[(i * W + j) * 3 + 2] = 255;
            }
			
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}