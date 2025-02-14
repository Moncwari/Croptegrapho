#include "Curve.h"
#include <ctime>
int main() {
    
    long long p = 13, a = 2, b = 2;
    clock_t start = clock();
    auto groups = findSimpleOrderGroups1(a, b, p);
    for (auto group : groups) {
        for (auto point : group) {
            std::cout << point << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Final";
    
    /*
    long long p = 974226079, a = 14, b = 13;
    long long x = -(p-1) / 2;
    long long x3 = powerMod(x, 3, p);
    long long y2 = (x3 + a * x + b) % p;
    std::cout << legendreSymbol(y2, p) << std::endl;
    long long y = sqrtMod(y2, p);
    Point q = Point(x , y, 14, p);
    std::cout << q << " " << BSGS(q) << std::endl;
    */
    

}