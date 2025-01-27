#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <set>
#include <cstdlib>
#include <ctime>


int findReverse(int x, int p) {
    for (int i = 0; i < p; ++i) {
        if ((i * x) % p == 1) return i; 
    }
    return -1;
}

int intPow(int base, int exp) {
    int result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result *= base;
        }
        base *= base;
        exp /= 2;
    }
    return result;
}

int legendreSymbol(int a, int p) {
    int ls = intPow(a, (p - 1) / 2) % p;
    return (ls == p - 1) ? -1 : ls;
}

int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

int lcm(int a, int b) {
    return (a / gcd(a, b)) * b;
}

long long findLCM(const std::vector<long long>& orders) {
    return std::accumulate(orders.begin(), orders.end(), 1LL, lcm);
}

std::vector<int> factorizeSlow(int n) {
    std::vector<int> factors;
    
    while (n % 2 == 0) {
        factors.push_back(2);
        n /= 2;
    }
    
    for (int i = 3; i * i <= n; i += 2) {
        while (n % i == 0) {
            factors.push_back(i);
            n /= i;
        }
    }

    if (n > 1) {
        factors.push_back(n);
    }
    
    return factors;
}

long long mulmod(long long a, long long b, long long mod) {
    long long res = 0;
    a %= mod;
    while (b > 0) {
        if (b % 2 == 1) res = (res + a) % mod;
        a = (a * 2) % mod;
        b /= 2;
    }
    return res;
}

long long f(long long x, long long c, long long n) {
    return (mulmod(x, x, n) + c) % n;
}

long long pollard_rho(long long n) {
    if (n % 2 == 0) return 2;
    long long x = 2, y = 2, c = 1;
    while (true) {
        x = f(x, c, n);
        y = f(f(y, c, n), c, n);
        long long d = std::__gcd(std::abs(x - y), n);
        if (d > 1) return d;
    }
}

std::vector<long long> factorize(long long n) {
    std::vector<long long> factors;
    while (n > 1) {
        if (n <= 1e9) {
            for (long long i = 2; i * i <= n; ++i) {
                while (n % i == 0) {
                    factors.push_back(i);
                    n /= i;
                }
            }
            if (n > 1) factors.push_back(n);
            break;
        }
        long long factor = pollard_rho(n);
        factors.push_back(factor);
        n /= factor;
    }
    return factors;
}

int powerMod(int a, int b, int p) {
    int res = 1;
    a %= p;
    while (b > 0) {
        if (b % 2 == 1) res = (1LL * res * a) % p;
        a = (1LL * a * a) % p;
        b /= 2;
    }
    return res;
}


long long sqrtMod_simple(long long a, long long p) {
    if (p % 4 != 3) return -1; 
    return powerMod(a, (p + 1) / 4, p);
}

long long sqrtMod_tonelli_shanks(long long a, long long p) {
    if (powerMod(a, (p - 1) / 2, p) != 1) return -1; 

    if (p % 4 == 3) return powerMod(a, (p + 1) / 4, p);

    long long q = p - 1, s = 0;
    while (q % 2 == 0) {
        q /= 2;
        ++s;
    }

    long long z = 2;
    while (powerMod(z, (p - 1) / 2, p) == 1) ++z;

    long long m = s;
    long long c = powerMod(z, q, p);
    long long t = powerMod(a, q, p);
    long long r = powerMod(a, (q + 1) / 2, p);

    while (t != 1) {
        long long i = 0, temp = t;
        while (temp != 1) {
            temp = (temp * temp) % p;
            ++i;
            if (i == m) return -1;
        }

        long long b = powerMod(c, (1LL << (m - i - 1)), p);
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
        m = i;
    }

    return r;
}

long long sqrtMod(long long a, long long p) {
    if (powerMod(a, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return powerMod(a, (p + 1) / 4, p);
    return sqrtMod_tonelli_shanks(a, p);
}

bool isQuadraticResidue(int n, int p) {
    if (n == 0) return true;
    int exp = (p - 1) / 2;
    return powerMod(n, exp, p) == 1;
}

struct Point {
    int x;
    int y;
    int a;
    int p;

    Point(int _x, int _y, int _a, int _p) {
        x = _x;
        y = _y;
        a = _a;
        p = _p;
    }

    void normalize(int &x, int p) {
        if (x > (p - 1) / 2) {
            x -= p;
        }
    }

    Point operator+(const Point &addend) {
        if (x == addend.x && y != addend.y) {
            return Point(0, 0, a, p);
        }
        int resX, resY;
        if (x != addend.x && y != addend.y) {
            resX = ((intPow((addend.y - y) * findReverse(addend.x - x, p), 2) - x - addend.x) % p + p) % p;
            normalize(resX, p);
            resY = ((((addend.y - y) * findReverse(addend.x - x, p)) * (x - resX) - y) % p + p) % p;
            normalize(resY, p);
        }
        if (x == addend.x && y == addend.y) {
            resX = ((intPow((3 * x * x + a) * findReverse(2 * y, p), 2) - 2 * x) % p + p) % p;
            normalize(resX, p);
            resY = ((((3 * x * x + a) * findReverse(2 * y, p)) * (x - resX) - y) % p + p) % p;
            normalize(resY, p);
        }
        return Point(resX, resY, a, p);
    }
    
    Point operator*(int scalar) {
        if (scalar == 0) return Point(0, 0, a, p);
        Point result(0, 0, a, p);
        Point base = *this;
        
        while (scalar > 0) {
            if (scalar % 2 == 1) {
                result = result + base;
            }
            base = base + base; 
            scalar /= 2;
        }
        return result;
    }

    bool operator==(const Point &other) const {
    return x == other.x && y == other.y && a == other.a && p == other.p;
    }

    friend std::ostream &operator<<(std::ostream &os, Point point) {
        os << '(' << point.x << ", " << point.y << ')';
        return os; 
    }

};

struct PointHash {
    std::size_t operator()(const Point& p) const {
        return std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1);
    }
};

std::map<int, int> fieldSquareRoots(int p) {
    std::vector<int> elems;
    std::map<int, int> squares;
    for (int i = 0; i <= p / 2; ++i) {
        elems.push_back(i);
    }
    for (int elem : elems) {
        int square = (elem * elem) % p;
        squares[square] = elem;
    }
    return squares;
}

std::vector<Point> buildCurve(int a, int b, int p) {
    std::vector<Point> points{Point(0, 0, a, p)};
    std::map<int, int> roots = fieldSquareRoots(p);
    for (int x = -p / 2; x <= p / 2; ++x) {
        int squareY = ((x * x * x + a * x + b) % p + p) % p;
        if (roots.find(squareY) == roots.end()) {
            continue;
        }
        std::pair<int, int> yRoots = {roots[squareY], -roots[squareY]};
        if (yRoots.first == yRoots.second) {
            points.push_back(Point(x, yRoots.first, a, p));
            continue;
        }
        points.push_back(Point(x, yRoots.first, a, p));
        points.push_back(Point(x, yRoots.second, a, p));
    }
    return points;
}

long long shuf(Point P) {

    if (P.x == 0 && P.y == 0) return 1;

    Point slow = P;
    Point fast = P + P;
    long long i = 1, j = 2;

    while (true) {
        
        slow = slow + P; i++;
        fast = fast + P; fast = fast + P; j += 2;

        if (slow.x == fast.x && slow.y == fast.y) {
            return j - i;
        }

        if ((slow.x == 0 && slow.y == 0) || (fast.x == 0 && fast.y == 0)) {
            return std::min(i, j);
        }
    }
}

int nativeFindOrder(int a, int b, int p) {
    int count = 1; 
    for (int x = 0; x < p; ++x) {
        int squareY = (x * x * x + a * x + b) % p;
        squareY = (squareY + p) % p;
        int legendre = legendreSymbol(squareY, p);
        if (legendre == 1) count += 2;
        else if (legendre == 0) count += 1;
    }
    return count;
}

int BSGS(Point P) {
    int p = P.p;
    if (P.x == 0 && P.y == 0) return 1;
    Point Q = P * (p + 1);
    long long m = static_cast<long long>(std::pow(p, 0.25) + 1);

    std::unordered_map<Point, long long,  PointHash> babySteps;

    Point temp = P;
    babySteps[Point(0, 0, P.a, p)] = 0;
    for (long long j = 1; j <= m; ++j) {
        babySteps[temp] = j;
        temp = temp + P;
    }

    Point giantStep = P * m;
    giantStep = giantStep + giantStep;
    long long l;
    long long M;
    for (long long k = - m; k <= m; ++k) {
        temp = Q + giantStep * k;
        Point revtemp = Point(temp.x, -temp.y, temp.a, temp.p);
        if (babySteps.find(temp) != babySteps.end()) {
            l = babySteps[temp];
            M = p + 1 + 2 * m * k;
            break;
        }
        if (babySteps.find(revtemp) != babySteps.end()) {
            l = -babySteps[revtemp];
            M = p + 1 + 2 * m * k;
            break;
        }
    }
    exit_loop:
    std::vector<long long> factors = factorize(M);
    for (long long f : factors) {
        temp = P * (M / f);
        if (temp.x == 0 && temp.y == 0) {
            M = M / f;
            goto exit_loop;
        }
    }

    return M; 
}

void removeInvalidOrders(long long lcmValue, std::set<long long>& orders) {
    for (auto it = orders.begin(); it != orders.end(); ) {
        if (lcmValue % *it != 0) {
            it = orders.erase(it);
        } else {
            ++it;
        }
    }
}

long long findOrder(long long a, long long b, long long p) {
    long long lowerBound = static_cast<long long>(std::ceil(p + 1 - 2 * std::sqrt(p)));
    long long upperBound = static_cast<long long>(std::floor(p + 1 + 2 * std::sqrt(p)));
    std::set<long long> possibleOrders; 
    std::cout << possibleOrders.size();
    for (long long N = lowerBound; N <= upperBound; ++N) {
        possibleOrders.insert(N);
    }
    std::vector<long long> pointOrders;
    std::srand(std::time(0));
    while (possibleOrders.size() != 1) {
        long long x = lowerBound + std::rand() % (upperBound - lowerBound + 1);
        long long y = sqrtMod(((x * x * x + a * x + b) % p + p) % p, p);
        if (y == -1) continue;
        Point P(x, y, a, b);
        std::cout << '(' << P.x << ", " << P.y << ')' << '\n';
        pointOrders.push_back(BSGS(P));
        std::cout << '-' << BSGS(P) << '-' << '\n';
        long long LCM = findLCM(pointOrders);
        removeInvalidOrders(LCM, possibleOrders);
    }
    return *possibleOrders.begin();
    
}

int main() {
    Point P(-2, -4, 6, 11);
    Point P1(-4, 5, 6, 11);
    
}