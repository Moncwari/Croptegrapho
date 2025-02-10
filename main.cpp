#include "Curve.hpp"

#define ll long long


int main() {
    ll p, a, b;
    std::cout << "Enter p: "; std::cin >> p;
    std::cout << "Enter a: "; std::cin >> a;
    std::cout << "Enter b: "; std::cin >> b;
    int mode;
    std::cout << "Enter 0 if you want to find order of curve,\n1 if you want buid curve,\n2 if you want find power of point,\n3 if you want find primary groups: \n"; std::cin >> mode;
    clock_t start = clock();
    if (mode == 0) {
        std::cout << "Order of curve: " << findOrder(a, b, p) << std::endl;
        std::cout << "Time: " << (clock() - start) / CLOCKS_PER_SEC << std::endl;
    } else if (mode == 1) {
        std::vector<Point> curve = buildCurve(a, b, p);
        int i = 0;
        for (auto point : curve) {
            i++;
            std::cout << point << std::endl;
            if (i == 50)
            {
                std::cout << "..." << std::endl;
                std::cout << curve.size() << " points at all" << std::endl;
                return 0;
            }
        }
        std::cout << curve.size() << " points at all" << std::endl;
    } else if (mode == 2) {
        std::cout << "Enter x: "; ll x; std::cin >> x;
        std::cout << "Enter y: "; ll y; std::cin >> y; normalize(y, p);
        ll y1 = sqrtMod(((intPow(x, 3, p) % p + a * x + b) % p + p) % p, p);
        normalize(y1, p);
        if (y != y1 && y != -y1) {
            std::cout << "No such point" << std::endl;
            return 0;
        }
        Point P(x, y, a, p);
        std::cout << "Enter power: "; ll k; std::cin >> k;
        Point Q = P * k;
        std::cout << "Your point: " << Q << std::endl;
    } else if (mode == 3) {
        std::vector<std::vector<Point>> groups = findSimpleOrderGroups(a, b, p);
        for (auto group : groups) {
            std::cout << "Size: " << group.size() << "\n";
            int i = 0;
            for (auto point : group) {
                std::cout << point << " ";
                ++i;
                if (i == 10)
                    break;
            }
            std::cout << std::endl;
        }
    }
}

