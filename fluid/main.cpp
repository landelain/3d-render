#include "save_svg.cpp"
#include <random>
#include "lbfgs.c"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <sstream>

#define PI 3.14159
#define VOL_FLUID 0.6

/*

    The basic outline of my code and most of the classes along with their methods 
    follow the implementation from the professor during the live-coding.

*/


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unif(0.0, 1.0);

class Voronoi_Diagram {
public : 
    
    std::vector<Polygon> diagram;
    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<Vector> unit_disk;
    int N_disk;

    Voronoi_Diagram(){
        N_disk = 20;
        unit_disk.resize(N_disk);
        for(int i = 0; i < N_disk; i++){
            double theta = i*(2.*PI / (double)N_disk);
            unit_disk[i] = Vector(-cos(theta), sin(theta), 0);
        }
    };

    Voronoi_Diagram(std::vector<Vector> points, int N_disk = 20){
        this->points = points;
        weights = std::vector<double>(points.size(), 0.1);
        
        this->N_disk = N_disk;
        unit_disk.resize(N_disk);
        for(int i = 0; i < N_disk; i++){
            double theta = i*(2.*PI / (double)N_disk);
            unit_disk[i] = Vector(-cos(theta), sin(theta), 0);
        }

    }

    Polygon Suth_Hod_clip(Vector& u, Vector& v, Polygon& Poly){
        Polygon result;
        Vector A,B;
        Vector N(u[1]- v[1], v[0] - u[0], 0);
        int n = Poly.vertices.size();

        for(int i = 0; i < n; i++){
            A = Poly.vertices[(i==0)? n-1: i-1];
            B = Poly.vertices[i];
            Vector P = A + dot((u-A), N)/dot((B-A), N) * (B-A);

            // B inside
            if(dot((u-B), N) >= 0){
                // A outside
                if(dot((u-A), N) < 0){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            // A inside
            else if (dot((u-A), N) >= 0){
                result.vertices.push_back(P);
            }
        }
        return result;
    }


    Polygon Voronoi_clip(Vector& Po, Vector& Pi, Polygon& Poly){
        Polygon result;
        Vector M = (Po + Pi)/2;
        Vector A,B;
        int n = Poly.vertices.size();

        for(int i = 0; i < n; i++){
            A = Poly.vertices[(i==0)? n-1: i-1];
            B = Poly.vertices[i];
            double t = dot((M-A), (Pi - Po))/dot((B-A), (Pi-Po));
            Vector P = A + t*(B-A);

            // B inside
            if((Po - B).norm2() <= (Pi - B).norm2()){
                // A outside
                if((Po - A).norm2() > (Pi - A).norm2()){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            // A inside
            else if ((Po - A).norm2() <= (Pi - A).norm2()){
                result.vertices.push_back(P);
            }
        }

        return result;
    }


    Polygon Laguerre_clip(Vector& Po, Vector& Pi, Polygon& Poly, double& Wo, double& Wi){
        Polygon result;
        Vector M = (Po + Pi)/2;
        Vector Mprime = M + (Wo - Wi)/(2. * (Po - Pi).norm2()) * (Pi - Po);
        Vector A,B;
        int n = Poly.vertices.size();

        for(int i = 0; i < n; i++){
            A = Poly.vertices[(i==0)? n-1: i-1];
            B = Poly.vertices[i];
            double t = dot((Mprime-A), (Pi - Po))/dot((B-A), (Pi-Po));
            Vector P = A + t * (B - A);
           
            // B inside
            if((Po - B).norm2() - Wo <= (Pi - B).norm2() - Wi + 1e-9){
                // A outside
                if((Po - A).norm2() - Wo > (Pi - A).norm2() - Wi + 1e-9){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            // A inside
            else if ((Po - A).norm2() - Wo <= (Pi - A).norm2() - Wi + 1e-9){
                result.vertices.push_back(P);
            }
        }

        return result;
    }


    void LLoyd(){
        std::vector<Vector> Centroids = {};

        for(int i = 0 ; i < diagram.size(); i++){
            Vector centroid = diagram[i].centroid();
            Centroids.push_back(centroid);
        }  
        this->points = Centroids;
    }

    void compute(){
        Polygon square;
        double size = 1;
        square.vertices.push_back(Vector(0, 0, 0));
        square.vertices.push_back(Vector(0, size, 0));
        square.vertices.push_back(Vector(size, size, 0));
        square.vertices.push_back(Vector(size, 0, 0));

        diagram.resize(points.size());
        if(weights.empty()){
            weights = std::vector<double>(points.size()+1, 1.);
        }

        #pragma omp parallel for schedule(dynamic, 1)
        for(int i = 0; i < points.size(); i++){
            Polygon V = square;
            for(int j = 0; j < points.size(); j++){
                if(i == j){continue;}

                V = Laguerre_clip(points[i], points[j], V, weights[i], weights[j]);
            }

            double radius = sqrt(weights[i] - weights[weights.size()-1]);
            //std::cout << radius << std::endl;

            for(int j = 0; j < N_disk; j++){
                int k = (j == N_disk - 1) ? 0 : j+1;
                Vector u = unit_disk[j]*radius + points[i];
                Vector v = unit_disk[k]*radius + points[i];
                V = Suth_Hod_clip(u, v, V);
            }

            diagram[i] = V;
        }
    }
};

static lbfgsfloatval_t evaluate(
     void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);

static int progress(
    void *instance, 
    const lbfgsfloatval_t *x, 
    const lbfgsfloatval_t *g, 
    const lbfgsfloatval_t fx, 
    const lbfgsfloatval_t xnorm, 
    const lbfgsfloatval_t gnorm, 
    const lbfgsfloatval_t step, 
    int n, 
    int k, 
    int ls
);

class OptimalTransport{
public: 

    Voronoi_Diagram diagram;

    OptimalTransport(){}

    OptimalTransport(Voronoi_Diagram diagram){
        this->diagram = diagram;
    }

    void optimize(){
        int i, ret = 0;
        lbfgsfloatval_t fx = 0.0;
        int N = diagram.weights.size();
        std::vector<lbfgsfloatval_t> weights(N);
        std::copy(diagram.weights.begin(), diagram.weights.end(), weights.begin());

        // std::cout << "maybe here " << std::endl;
        // std::cout << diagram.weights[100] << std::endl;

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        param.epsilon = 1e-8;
        param.max_iterations = 1000;

        ret = lbfgs(N, &weights[0], &fx, evaluate, NULL, (void*) this, &param);

        std::copy(weights.begin(), weights.begin() + (N - 1), diagram.weights.begin());
        //diagram.weights = weights;
        //memcpy(&diagram.weights[0]);

        diagram.compute();
        
    }

};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
){

    OptimalTransport* ot = static_cast<OptimalTransport*>(instance);
    lbfgsfloatval_t fx = 0.0;
    int N = ot->diagram.weights.size();
    memcpy(&ot->diagram.weights[0], x, N*sizeof(x[0]));

    //ot->diagram.weights.assign(x, x + n); //memcpy
    //memcpy(&ot->diagram.weights[0], x, n*sizeof(x[0]));
    ot->diagram.compute();
    double sum_fluid = 0;

    for(int i=0; i < N-1; i++){
        double current_area = ot->diagram.diagram[i].area();
        sum_fluid += current_area;
        double area = VOL_FLUID/(N-1);

        g[i] = -( area - current_area);
        fx += ot->diagram.diagram[i].integral_sq_dist(ot->diagram.points[i]) - x[i] * (current_area - area) ;
    }

    double estimated_air_vol = 1 - sum_fluid;
    double desired_air_vol = (1. - VOL_FLUID);
    g[N-1] = -(desired_air_vol - estimated_air_vol);
    fx+= x[N-1] * (desired_air_vol - estimated_air_vol);

    return -fx;

}

static int progress(
    void *instance, 
    const lbfgsfloatval_t *x, 
    const lbfgsfloatval_t *g, 
    const lbfgsfloatval_t fx, 
    const lbfgsfloatval_t xnorm, 
    const lbfgsfloatval_t gnorm, 
    const lbfgsfloatval_t step, 
    int n, 
    int k, 
    int ls
){

    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;

}


int sgn(double x){
    if( x> 0){return 1;}
    if(x < 0){return -1;}
    return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    //if (i < N) {   // the N first particles may represent fluid, displayed in blue
                    //	image[((H - y - 1)*W + x) * 3] = 0;
                    //	image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    //	image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    //}
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
};

class Fluid {
public:

    OptimalTransport ot;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    int N;

    Fluid(int N){
        this->N = N;
        particles.resize(N);
        velocities.resize(N, Vector(0, 0, 0));
        for(int i = 0; i < N; i++){
            particles[i] = Vector(unif(gen), unif(gen), 0);
        }
        this->ot.diagram.points = particles;
        this->ot.diagram.weights.resize(N+1);
        std::fill(this->ot.diagram.weights.begin(), this->ot.diagram.weights.end(), 0.);
        this->ot.diagram.weights[N] = 0.9;

    }

    Fluid(int N, OptimalTransport ot){
        this->N = N;
        this->ot = ot;

        particles.resize(N);
        velocities.resize(N, Vector(0, 0, 0));
        for(int i = 0; i < N; i++){
            particles[i] = Vector(unif(gen), unif(gen), 0);
        }
        
        ot.diagram.points = particles;
        ot.diagram.weights.resize(N+1);
        std::fill(ot.diagram.weights.begin(), ot.diagram.weights.end(), 1.);

    }   

    void time_step(double dt){

        ot.optimize();

        // std::cout << "here look" << std::endl;
        // std::cout << ot.diagram.weights.size() << std::endl;
        // std::cout << ot.diagram.weights[100] << std::endl;

        double eps = 0.004;
        double eps2 = eps*eps;
        Vector g(0, -9.81, 0);
        double mass = 200;

        for(int i = 0; i < particles.size(); i++){

            Vector center_cell = ot.diagram.diagram[i].centroid() ;// define centroid
            Vector spring_force = (center_cell - particles[i])/(eps2);
            Vector all_forces = g*mass + spring_force;

            velocities[i] = velocities[i] +  dt/mass * all_forces;
            //std::cout << particles[i][1] << std::endl;
            particles[i] = particles[i] + dt * velocities[i];
            //std::cout << particles[i][1] << std::endl;
        }

        ot.diagram.points = particles;

    }

    void simulation(int pics= 100){
        double dt = 0.005;
        for(int i = 0; i < pics ; i ++){
            std::cout << i+1 << " over " << pics << std::endl;
            time_step(dt);
            save_frame(ot.diagram.diagram, "test2/animation", i);
        }
    }

};

int main(){

    Fluid fluid(100);
    fluid.simulation(100);
    
    // int n_points = 100;
    
    // std::vector<Vector> points;
    // for(int i = 0; i < n_points ; i++){
    //     double x = unif(gen);
    //     double y = unif(gen);
    //     Vector point(x, y, 0);
    //     points.push_back(point);
    // }

    // Voronoi_Diagram diagram(points);
    // OptimalTransport ot(diagram);
    // ot.diagram.compute();
    // ot.optimize();
    // //std::cout << "hey" << std::endl;
    // //ot.diagram.compute();
    // //std::cout << "hey" << std::endl;

    // save_svg(ot.diagram.diagram, "firstp.svg", "lightblue");
    //std::cout << "hey" << std::endl;

    // for(int i = 0; i< n_points ; i++){
    //     std::cout << ot.diagram.diagram[i].vertices[0][0] << ", " << ot.diagram.diagram[i].vertices[0][1] << std::endl;
    // }

    // for(int i = 0; i< 20 ; i++){
    //     std::cout << diagram.points[i][0] << ", " << diagram.points[i][0] << std::endl;
    // }

    // for(int i = 0; i < 10; i++){
    //     ot.diagram.LLoyd();
    //     //std::cout << "hey" << std::endl;
    //     ot.diagram.compute();
    //     //std::cout << "hey" << std::endl;
    // }
    
    
    // std::cout << diagram.diagram.size() << std::endl;
    // std::cout << diagram.points.size() << std::endl;
    //std::cout << " " << std::endl;

    // for(int i = 0; i< 20 ; i++){
    //     std::cout << diagram.diagram[i].vertices[0][0] << ", " << diagram.diagram[i].vertices[0][1] << std::endl;
    // }

    //save_svg(ot.diagram.diagram, "secondp.svg", "lightblue");

    return 0;
}