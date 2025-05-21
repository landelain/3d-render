#include "save_svg.cpp"
#include <random>
#include "lbfgs.c"


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unif(0.0, 1.0);

class Voronoi_Diagram {
public : 
    
    std::vector<Polygon> diagram;
    std::vector<Vector> points;
    std::vector<double> weights;

    Voronoi_Diagram(){};

    Voronoi_Diagram(std::vector<Vector> points){
        this->points = points;
        weights = std::vector<double>(points.size(), 0.1);
        // for(int i = 0; i < points.size(); i++){
        //     weights[i] = unif(gen)*0.01;
        // }
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
            if((Po - B).norm2() - Wo <= (Pi - B).norm2() - Wi){
                // A outside
                if((Po - A).norm2() - Wo > (Pi - A).norm2() - Wi){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            // A inside
            else if ((Po - A).norm2() - Wo <= (Pi - A).norm2() - Wi){
                result.vertices.push_back(P);
            }
        }

        return result;
    }


    void LLoyd(){
        std::vector<Vector> Centroids = {};

        for(int i = 0 ; i < diagram.size(); i++){
            double sum = 0;
            int n = diagram[i].vertices.size();

            for(int j = 0 ; j < n; j++){
                sum += (diagram[i].vertices[j][0] * diagram[i].vertices[(j == n-1) ? 0 : j+1][1]) - (diagram[i].vertices[j][1] * diagram[i].vertices[(j == n-1) ? 0 : j+1 ][0]);
            }
            double A = 1./2. * sum;
            sum = 0;

            for(int j = 0 ; j < n; j++){
                sum += (diagram[i].vertices[j][0] + diagram[i].vertices[(j == n-1) ? 0 : j+1][0]) * (diagram[i].vertices[j][0] * diagram[i].vertices[(j == n-1) ? 0 : j+1][1] - diagram[i].vertices[j][1] * diagram[i].vertices[(j == n-1) ? 0 : j+1 ][0]);
            }
            double Cx = 1./(6. * A) * sum;
            sum = 0;

            for(int j = 0 ; j < n; j++){
                sum += (diagram[i].vertices[j][1] + diagram[i].vertices[(j == n-1) ? 0 : j+1][1]) * ( diagram[i].vertices[j][0] * diagram[i].vertices[(j == n-1) ? 0 : j+1][1] - diagram[i].vertices[j][1] * diagram[i].vertices[(j == n-1) ? 0 : j+1 ][0]);
            }
            double Cy = 1./(6. * A) * sum;

            Centroids.push_back(Vector(Cx, Cy, 0.));
        }  

        this->points = Centroids;
    }


    void compute(){
        diagram = {};
        Polygon square;
        double size = 1;
        square.vertices.push_back(Vector(0, 0, 0));
        square.vertices.push_back(Vector(size, 0, 0));
        square.vertices.push_back(Vector(size, size, 0));
        square.vertices.push_back(Vector(0, size, 0));

        for(int i = 0; i < points.size(); i++){
            Polygon V = square;
            for(int j = 0; j < points.size(); j++){
                if(i == j){continue;}
                //V = Voronoi_clip(points[i], points[j], V);
                V = Laguerre_clip(points[i], points[j], V, weights[i], weights[j]);
            }

            diagram.push_back(V);
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
        int N = diagram.points.size();
        std::vector<lbfgsfloatval_t> weights(N, 0.0);

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        param.epsilon = 1e-8;
        param.max_iterations = 1000;

        ret = lbfgs(N, &weights[0], &fx, evaluate, progress, (void*) this, &param);

        diagram.weights = weights;

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

    ot->diagram.weights.assign(x, x + n);
    ot->diagram.compute();
    
    for(int i=0; i < n; i++){
        double current_area = ot->diagram.diagram[i].area();
        double area = 1./n;

        g[i] = -(area - current_area);
        fx += ot->diagram.diagram[i].integral_sq_dist(ot->diagram.points[i]) - x[i] * (current_area - area) ;
    }

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

int main(){

    int n_points = 100;
    
    std::vector<Vector> points;
    for(int i = 0; i < n_points ; i++){
        double x = unif(gen);
        double y = unif(gen);
        Vector point(x, y, 0);
        points.push_back(point);
    }

    Voronoi_Diagram diagram(points);
    OptimalTransport ot(diagram);
    ot.diagram.compute();
    ot.optimize();
    //std::cout << "hey" << std::endl;
    //ot.diagram.compute();
    //std::cout << "hey" << std::endl;

    save_svg(ot.diagram.diagram, "firstp.svg", "lightblue");
    //std::cout << "hey" << std::endl;

    // for(int i = 0; i< n_points ; i++){
    //     std::cout << ot.diagram.diagram[i].vertices[0][0] << ", " << ot.diagram.diagram[i].vertices[0][1] << std::endl;
    // }

    // for(int i = 0; i< 20 ; i++){
    //     std::cout << diagram.points[i][0] << ", " << diagram.points[i][0] << std::endl;
    // }

    for(int i = 0; i < 10; i++){
        ot.diagram.LLoyd();
        //std::cout << "hey" << std::endl;
        ot.diagram.compute();
        //std::cout << "hey" << std::endl;
    }
    
    
    // std::cout << diagram.diagram.size() << std::endl;
    // std::cout << diagram.points.size() << std::endl;
    //std::cout << " " << std::endl;

    // for(int i = 0; i< 20 ; i++){
    //     std::cout << diagram.diagram[i].vertices[0][0] << ", " << diagram.diagram[i].vertices[0][1] << std::endl;
    // }

    save_svg(ot.diagram.diagram, "secondp.svg", "lightblue");

    return 0;
}