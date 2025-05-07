#include "save_svg.cpp"
#include <random>


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unif(0.0, 1.0);

class Voronoi_Diagram {
public : 
    
    std::vector<Polygon> diagram;
    std::vector<Vector> points;

    Voronoi_Diagram(){};

    Voronoi_Diagram(std::vector<Vector> points){
        this->points = points;
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
    
    void compute(){
        Polygon square;
        square.vertices.push_back(Vector(0, 0, 0));
        square.vertices.push_back(Vector(1, 0, 0));
        square.vertices.push_back(Vector(1, 1, 0));
        square.vertices.push_back(Vector(0, 1, 0));

        for(int i = 0; i < points.size(); i++){
            Polygon V = square;
            for(int j = 0; j < points.size(); j++){
                if(i == j){continue;}
                V = Voronoi_clip(points[i], points[j], V);
            }

            diagram.push_back(V);

        }

    }
};

int main(){
    
    std::vector<Vector> points;
    for(int i = 0; i < 30 ; i++){
        double x = unif(gen);
        double y = unif(gen);
        Vector point(x, y, 0);
        points.push_back(point);
    }

    Voronoi_Diagram diagram(points);
    diagram.compute();
    save_svg(diagram.diagram, "test.svg");

    return 0;
}