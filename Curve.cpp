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
#include <random>


int findReverse(int x, int p) {
    if (x < 0) x = (x % p + p) % p;
    int a = x, b = p, u = 1, v = 0;
    while (b != 0) {
        int t = a / b;
        a -= t * b; std::swap(a, b);
        u -= t * v; std::swap(u, v);
    }
    if (a != 1) return -1;
    return (u % p + p) % p;
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

long long powerMod(long long a, long long b, long long p) {
    long long res = 1;
    a %= p;
    while (b) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

bool isPrimeFermat(long long n, int iterations = 5) {
    if (n < 4)
        return (n == 2 || n == 3); 

    srand(time(0));

    for (int i = 0; i < iterations; i++) {
        long long a = 2 + rand() % (n - 3);
        if (powerMod(a, n - 1, n) != 1)
            return false;
    }
    return true;
}

bool millerTest(long long d, long long n) {
    if (n < 4) return (n == 2 || n == 3);
    
    long long a = 2 + rand() % (n - 4);
    long long x = powerMod(a, d, n);

    if (x == 1 || x == n - 1)
        return true;

    while (d != n - 1) {
        x = (__int128_t(x) * x) % n;
        d *= 2;

        if (x == 1) return false;
        if (x == n - 1) return true;
    }

    return false;
}

bool isPrimeMillerRabin(long long n, int iterations = 10) {
    if (n <= 4) return (n == 2 || n == 3);

    long long d = n - 1;
    while (d % 2 == 0)
        d /= 2;

    for (int i = 0; i < iterations; i++) {
        if (!millerTest(d, n))
            return false;
    }
    return true;
}

long long f(long long x, long long c, long long n) {
    return (powerMod(x, x, n) + c) % n;
}

long long pollardRho(long long n) {
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
        long long factor = pollardRho(n);
        factors.push_back(factor);
        n /= factor;
    }
    return factors;
}

long long sqrtModTonnelliShanks(long long a, long long p) {
    if (a == 0) return 0;
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
    if (a == 0) return 0; 
    if (powerMod(a, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return powerMod(a, (p + 1) / 4, p);
    return sqrtModTonnelliShanks(a, p);
}

struct Point {
    int x, y, a, p, order;

    Point() : x(0), y(0), a(0), p(-1) {}

    Point(int _x, int _y, int _a, int _p) : x(_x), y(_y), a(_a), p(_p) {}

    static void normalize(int &x, int p) {
        
        if (x > (p - 1) / 2) {
            x -= p;
        }
    }

    Point operator+(const Point &addend) {
        if (x == addend.x && y != addend.y) {
            return Point();
        }
        if (p == -1) return addend;
        if (addend.p == -1) return *this;
        int resX, resY;
        if (x == addend.x && y == addend.y) {
            resX = ((intPow((3 * x * x + a) * findReverse(2 * y, p), 2) - 2 * x) % p + p) % p;
            normalize(resX, p);
            resY = ((((3 * x * x + a) * findReverse(2 * y, p)) * (x - resX) - y) % p + p) % p;
            normalize(resY, p);
        }
        else {
            resX = ((intPow((addend.y - y) * findReverse(addend.x - x, p), 2) - x - addend.x) % p + p) % p;
            normalize(resX, p);
            resY = ((((addend.y - y) * findReverse(addend.x - x, p)) * (x - resX) - y) % p + p) % p;
            normalize(resY, p);
        }
        return Point(resX, resY, a, p);
    }

    Point operator*(int scalar) const {
        if (scalar == 0 || p == -1) return Point();
        Point result = *this;
        Point base = *this;
        if (scalar == -1){
            return Point(x, -y, a, p);
        }
        else if (scalar < 0){
            return  result * (-scalar) * (-1);
        }
        scalar--;
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

    bool operator<(const Point &other) const {
    return std::tie(x, y, a, p) < std::tie(other.x, other.y, a, p);
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &point) {
        if (point.p == -1) os << '0';
        else os << '(' << point.x << ", " << point.y << ')';
        return os;
    }
};

class Polynomial {

  std::vector<int> coefficients;
                                 
  int galuaDet;

public:

  void setDet(int det) {
    galuaDet = det;
  }

  int getDet() const { return galuaDet; }

  int getCoef(const int degree) const { return coefficients[degree]; }

  void trim() {
    while (coefficients.size() > 1 && coefficients.back() == 0) {
      coefficients.pop_back();
    }
  }

  void modulo() {
    for (size_t i = 0; i < coefficients.size(); ++i) {
      coefficients[i] = (coefficients[i] % galuaDet + galuaDet) % galuaDet;
    }
    trim();
  }

  Polynomial() {
    coefficients = std::vector<int>({0});
    galuaDet = INT_MAX;
  }

  Polynomial(const std::vector<int> &coeffs, int det = INT_MAX)
      : galuaDet{det},
        coefficients{coeffs.rbegin(), std::reverse_iterator{std::find_if(
                                          coeffs.begin(), coeffs.end() - 1,
                                          [](int a) { return a; })}} {}
  int degree() const { return coefficients.size() - 1; }

  const int galua_division(const int divident, const int divider) const {
    if (divider == 1)
      return divident;

    if (galuaDet != INT_MAX) {
      for (int i = 1; i < galuaDet; ++i) {
        if (i * divider % galuaDet == 1) {
          return (divident * i % galuaDet);
        }
      }
    }
    
    return -1;
  }

  Polynomial operator+(const Polynomial &addend) const {
    size_t sumDegree =
        std::max(coefficients.size(), addend.coefficients.size());
    std::vector<int> sum(sumDegree, 0);

    for (size_t i = 0; i < sumDegree; ++i) {
      int a = (i < coefficients.size()) ? coefficients[i]
                                        : 0;
      int b = (i < addend.coefficients.size()) ? addend.coefficients[i] : 0;
      sum[i] = a + b;
    }
    std::reverse(sum.begin(), sum.end());
    Polynomial res = Polynomial(sum, galuaDet);
    res.modulo();
    res.trim();
    return res;
  }

  Polynomial operator-(const Polynomial &subtrahend) const {
    size_t diffDegree =
        std::max(coefficients.size(), subtrahend.coefficients.size());
    std::vector<int> difference(diffDegree, 0);

    for (size_t i = 0; i < diffDegree; ++i) {
      int a = (i < coefficients.size()) ? coefficients[i]
                                        : 0;
      int b = (i < subtrahend.coefficients.size())
                  ? subtrahend.coefficients[i]
                  : 0;
      difference[i] = a - b;
    }
    std::reverse(difference.begin(), difference.end());
    Polynomial ans(difference, galuaDet);
    ans.modulo();
    ans.trim();
    return ans;
  }

  Polynomial operator*(const Polynomial &cofactor) const {
    size_t compDegree = coefficients.size() + cofactor.coefficients.size() - 1;
    std::vector<int> comp(compDegree, 0);

    for (size_t i = 0; i < coefficients.size(); ++i) {
      for (size_t j = 0; j < cofactor.coefficients.size(); ++j) {
        comp[i + j] +=
            coefficients[coefficients.size() - i - 1] *
            cofactor.coefficients[cofactor.coefficients.size() - j - 1];
      }
    }

    return Polynomial(comp, galuaDet);
  }

  Polynomial operator%(const Polynomial &polydivisor) const {
    std::vector<int> divident(coefficients.begin(), coefficients.end());
    std::vector<int> divisor(polydivisor.coefficients.begin(),
                             polydivisor.coefficients.end());
    std::reverse(divident.begin(), divident.end());
    std::reverse(divisor.begin(), divisor.end());
    std::vector<int> remainder = divident;
    if (divisor.empty() || (divisor.size() == 1 && divisor[0] == 0)) {
      throw std::invalid_argument("Division by zero polynomial");
    }
    int remainderSize = remainder.size();
    int divisorSize = divisor.size();
    int k = 0;
    while (remainderSize >= divisorSize) {
      int coeff = galua_division(remainder[k], divisor[0]);
      for (size_t i = 0; i < divisor.size(); ++i) {
        remainder[k + i] -= coeff * divisor[i];
      }
      k += 1;
      remainderSize -= 1;
    }

    std::vector<int> intremainder(remainder.end() - remainderSize,
                                  remainder.end());
    Polynomial c(intremainder, polydivisor.getDet());
    c.modulo();
    return c;
  }

  bool operator==(const Polynomial &other) const {
    return (coefficients == other.coefficients && galuaDet == other.galuaDet);
  }

  friend std::ostream &operator<<(std::ostream &os, Polynomial poly1) {
    const std::string degreeSymbols[] = {"",  "",  "²", "³", "⁴",
                                         "⁵", "⁶", "⁷", "⁸", "⁹"};
    Polynomial poly = poly1;
    std::reverse(poly.coefficients.begin(), poly.coefficients.end());
    size_t degree = poly.coefficients.size() - 1;
    bool isFirst = true;

    for (size_t i = 0; i < poly.coefficients.size(); ++i) {
      int coeff = poly.coefficients[i];
      size_t currentDegree = degree - i;

      if (coeff != 0) {
        if (!isFirst) {
          os << (coeff > 0 ? " + " : " - ");
        } else if (coeff < 0) {
          os << "-";
        }

        if (std::abs(coeff) != 1 || currentDegree == 0) {
          os << std::abs(coeff);
        }

        if (currentDegree > 0) {
          os << "x";
          if (currentDegree < 10) {
            os << degreeSymbols[currentDegree];
          } else {
            os << "^" << currentDegree;
          }
        }

        isFirst = false;
      }
    }

    if (isFirst) {
      os << "0";
    }

    return os;
  }
};

namespace std {
    template <>
    struct hash<Point> {
        size_t operator()(const Point &point) const {
            size_t h1 = std::hash<int>()(point.x);
            size_t h2 = std::hash<int>()(point.y);
            return h1 ^ (h2 * 19); 
        }
    };
}

static void normalize(int& coord, int p) {
    if (coord > (p - 1) / 2) {
        coord -= p;
    }
    if (coord < -(p - 1) / 2) {
        coord += p;
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
    Point Q = P * (p + 1);
    if (P == Point()){
        return 1;
    } 

    long long m = static_cast<long long>(std::pow(p, 0.25) + 1);

    std::unordered_map<Point, long long> babySteps;
    Point temp = P;
    for (long long j = 1; j <= m; ++j) {
        if (temp == Point()){
            return j;
        }
        babySteps[temp] = j;
        temp = temp + P;
    }
    Point giantStep = P * m;
    long long l;
    long long M;
    for (long long k = - 2 * m; k <= 2 * m; ++k) {
        temp = Q + giantStep * k;
        Point revtemp = temp * (-1);
        if (babySteps.find(temp) != babySteps.end()) {
            l = babySteps[temp];   
            if (p + 1 + m * k - l != 0){
                M = p + 1 + m * k - l;
                break;
            }
        }
        if (babySteps.find(revtemp) != babySteps.end()) {
            l = babySteps[revtemp];
            if (p + 1 + m * k + l != 0){
                M = p + 1 + m * k + l;
                break;
            }
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

std::vector<Point> buildCurve(int a, int b, int p) {
    std::vector<Point> points{{0, 0, 0, -1}}; 
    for (int x = -p / 2; x <= p / 2; ++x) {
        int squareY = ((x * x * x + a * x + b) % p + p) % p;
        int s = sqrtMod(squareY, p);
        if (s == -1) {
            continue;
        }
        std::pair<int, int> yRoots{s, -s};
        if (yRoots.first == yRoots.second) {
            normalize(yRoots.first, p);
            points.emplace_back(x, yRoots.first, a, p);
            continue;
        }
        normalize(yRoots.first, p);
        normalize(yRoots.second, p);
        points.emplace_back(x, yRoots.first, a, p);
        points.emplace_back(x, yRoots.second, a, p);
    }

    return points;
}


void removeInvalidOrders(long long lcmValue, std::set<long long>& orders) {
    for (auto it = orders.begin(); it != orders.end(); ) {
        if (*it % lcmValue != 0) {
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
    
    for (long long N = lowerBound; N <= upperBound; ++N) {
        possibleOrders.insert(N);
    }

    std::vector<long long> pointOrders;
    std::srand(std::time(0));
    while (possibleOrders.size() > 1) {
        long long x = -p / 2 + std::rand() % p;
        long long y = sqrtMod(((x * x * x + a * x + b) % p + p) % p, p);
        if  (y > p/2) y -= p;
        if (y == -1) continue;
        Point P(x, y, a, p);
        pointOrders.push_back(BSGS(P));
        long long LCM = findLCM(pointOrders);
        removeInvalidOrders(LCM, possibleOrders);
    }
    return *possibleOrders.begin();
    
}

std::vector<std::vector<Point>> findSimpleOrderGroups(long long a, long long b, long long p) {
    std::vector<Point> curve = buildCurve(a, b, p);
    std::vector<std::vector<Point>> simpleOrderGroups;
    std::set<Point> gens;
    for (auto point : curve) {
        bool fl = false;
        long long n = BSGS(point);
        if (!isPrimeMillerRabin(n)) continue;
        std::vector<Point> group;
        Point cur{0, 0, 0, -1};
        for (long long i = 0; i < n; ++i) {
            group.push_back(cur);
            cur = cur + point;
            if (gens.find(cur) != gens.end())
            {
                fl = true;
                break;
            }
        }
        if (fl) continue;
        gens.insert(point);
        simpleOrderGroups.push_back(group);
    }
    return simpleOrderGroups;
}

Point randPoint(long long a, long long b, long long p) {
    while (true) {
        long long x = -p / 2 + std::rand() % p;
        long long y = sqrtMod(((x * x * x + a * x + b) % p + p) % p, p);
        if  (y > p/2) y -= p;
        if (y == -1) continue;
        Point P(x, y, a, p);
        return P;
    }
}

std::vector<long long> findPrimes(long long q) {
    long long pow = 4 * (sqrt(q) + 1);
    std::vector<long long> primes;
    long long n = 0;
    while (std::accumulate(primes.begin(), primes.end(), 1LL, std::multiplies<>()) < pow) {
        for (long long i = n + 1; true; ++i) {
            if (isPrimeFermat(i)) {
                primes.push_back(i);
                n = i;
                break;
            }
        }
    }
    return primes;
}

Polynomial shufPoly(long long l) {
    if (l == 0) return Polynomial();
    if (l == 1) return Polynomial({1});
    if (l == 2) return Polynomial({1, 0});
    else return Polynomial({1, 0}) * shufPoly(l - 1) - shufPoly(l - 2);
}

long long frobenius(long long l, long long p, long long a, long long b) {
    Polynomial poly = shufPoly(l);
    Point P = randPoint(a, b, p);
    
}

int main() {

}