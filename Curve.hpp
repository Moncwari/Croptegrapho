#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

long long powerMod(long long a, long long b, long long p);

long long findReverse(long long x, long long p) {
  if (x < 0) x = (x % p + p) % p;
  long long a = x, b = p, u = 1, v = 0;
  while (b != 0) {
      long long t = a / b;
      a -= t * b; std::swap(a, b);
      u -= t * v; std::swap(u, v);
  }
  if (a != 1) return -1;
  return (u % p + p) % p;
}
long long intPow(long long base, long long exp, long long p) {
  long long result = 1;
  base %= p;
  while (exp > 0) {
    if (exp % 2 == 1) {
      result *= base;
      result %= p;
    }
    base *= base;
    base %= p;
    exp /= 2;
    
  }
  return result;
}

long long legendreSymbol(long long a, long long p) {
  long long ls = intPow(a, (p - 1) / 2, p);
  return (ls == p - 1) ? -1 : ls;
}

long long gcd(long long a, long long b) {
  if (b == 0)
    return a;

  return gcd(b, a % b);
}

long long lcm(long long a, long long b) { return (a / gcd(a, b)) * b; }

long long findLCM(const std::vector<long long> &orders) {
  return std::accumulate(orders.begin(), orders.end(), 1LL, lcm);
}

std::vector<long long> factorizeSlow(long long n) {
  std::vector<long long> factors;

  while (n % 2 == 0) {
    factors.push_back(2);
    n /= 2;
  }

  for (long long i = 3; i * i <= n; i += 2) {
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
    if (b & 1)
      res = res * a % p;
    a = a * a % p;
    b >>= 1;
  }
  return res;
}

bool isPrimeFermat(long long n, long long iterations = 5) {
  if (n < 4)
    return (n == 2 || n == 3);

  srand(time(0));

  for (long long i = 0; i < iterations; i++) {
    long long a = 2 + rand() % (n - 3);
    if (powerMod(a, n - 1, n) != 1)
      return false;
  }
  return true;
}

bool millerTest(long long d, long long n) {
  if (n < 4)
    return (n == 2 || n == 3);

  long long a = 2 + rand() % (n - 4);
  long long x = powerMod(a, d, n);

  if (x == 1 || x == n - 1)
    return true;

  while (d != n - 1) {
    x = (__int128_t(x) * x) % n;
    d *= 2;

    if (x == 1)
      return false;
    if (x == n - 1)
      return true;
  }

  return false;
}

bool isPrimeMillerRabin(long long n, long long iterations = 10) {
  if (n <= 4)
    return (n == 2 || n == 3);

  long long d = n - 1;
  while (d % 2 == 0)
    d /= 2;

  for (long long i = 0; i < iterations; i++) {
    if (!millerTest(d, n))
      return false;
  }
  return true;
}

long long f(long long x, long long c, long long n) {
  return (powerMod(x, x, n) + c) % n;
}

long long pollardRho(long long n) {
  if (n % 2 == 0)
    return 2;
  long long x = 2, y = 2, c = 1;
  while (true) {
    x = f(x, c, n);
    y = f(f(y, c, n), c, n);
    long long d = std::__gcd(std::abs(x - y), n);
    if (d > 1)
      return d;
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
      if (n > 1)
        factors.push_back(n);
      break;
    }
    long long factor = pollardRho(n);
    factors.push_back(factor);
    n /= factor;
  }
  return factors;
}

long long sqrtModTonnelliShanks(long long a, long long p) {
  if (a == 0)
    return 0;
  if (powerMod(a, (p - 1) / 2, p) != 1)
    return -1;

  if (p % 4 == 3)
    return powerMod(a, (p + 1) / 4, p);

  long long q = p - 1, s = 0;
  while ((q & 1) == 0) {
    q >>= 1;
    ++s;
  }

  long long z = 2;
  while (powerMod(z, (p - 1) / 2, p) == 1)
    ++z;

  long long m = s;
  long long c = powerMod(z, q, p);
  long long t = powerMod(a, q, p);
  long long r = powerMod(a, (q + 1) / 2, p);

  while (t != 1) {
    long long i = 0, temp = t;
    for (; i < m && temp != 1; ++i)
      temp = (temp * temp) % p;

    if (i == m)
      return -1;

    long long b = powerMod(c, 1LL << (m - i - 1), p);
    c = (b * b) % p;
    t = (t * c) % p;
    r = (r * b) % p;
    m = i;
  }

  return r;
}

long long sqrtMod(long long a, long long p) {
  if (a == 0)
    return 0;
  if (powerMod(a, (p - 1) / 2, p) != 1)
    return -1;
  if (p % 4 == 3)
    return powerMod(a, (p + 1) / 4, p);
  return sqrtModTonnelliShanks(a, p);
}

struct Point {
  long long x, y, a, p, order;

  Point() : x(0), y(0), a(0), p(-1) {}

  Point(long long _x, long long _y, long long _a, long long _p)
      : x(_x), y(_y), a(_a), p(_p) {}

  static void normalize(long long &x, long long p) {
    x %= p;
    x += p;
    x %= p;
    if (x > (p - 1) / 2) {
      x -= p;
    }
  }

  Point operator+(const Point &addend) {
    if (x == addend.x && y != addend.y) {
      return Point();
    }
    if (p == -1)
      return addend;
    if (addend.p == -1)
      return *this;
    long long resX, resY;
    if (x == addend.x && y == addend.y) {
      resX =
          ((intPow((3 * x * x + a) * findReverse(2 * y, p), 2, p) - 2 * x) % p + p) % p;
      normalize(resX, p);
      resY = ((((3 * x * x + a) * findReverse(2 * y, p)) * (x - resX) - y) % p +
              p) %
             p;
      normalize(resY, p);
    } else {
      resX = ((intPow((addend.y - y) * findReverse(addend.x - x, p), 2, p) - x -
               addend.x) %
                  p +
              p) %
             p;
      normalize(resX, p);
      resY =
          ((((addend.y - y) * findReverse(addend.x - x, p)) * (x - resX) - y) %
               p +
           p) %
          p;
      normalize(resY, p);
    }
    return Point(resX, resY, a, p);
  }

  Point operator*(long long scalar) const {
    if (scalar == 0 || p == -1)
      return Point();
    Point result = *this;
    Point base = *this;
    if (scalar == -1) {
      return Point(x, -y, a, p);
    } else if (scalar < 0) {
      return result * (-scalar) * (-1);
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
    if (point.p == -1)
      os << '0';
    else
      os << '(' << point.x << ", " << point.y << ')';
    return os;
  }
};

namespace std {
template <> struct hash<Point> {
  size_t operator()(const Point &point) const {
    size_t h1 = std::hash<long long>()(point.x);
    size_t h2 = std::hash<long long>()(point.y);
    return h1 ^ (h2 * 19);
  }
};
} // namespace std

static void normalize(long long &coord, long long p) {
  if (coord > (p - 1) / 2) {
    coord -= p;
  }
  if (coord < -(p - 1) / 2) {
    coord += p;
  }
}

/**
 * @brief Computes the number of points on an elliptic curve over a finite
 * field.
 *
 * This function calculates the number of solutions (x, y) to the elliptic curve
 * equation y^2 = x^3 + ax + b over the finite field with prime order p.
 *
 * @param a The coefficient "a" of the elliptic curve equation.
 * @param b The coefficient "b" of the elliptic curve equation.
 * @param p The prime modulus of the elliptic curve.
 * @return The number of points on the elliptic curve.
 */
long long nativeFindOrder(long long a, long long b, long long p) {
  long long count = 1; // Start with the neutral element (point at infinity)

  // Iterate over all possible x values in the field F_p
  for (long long x = 0; x < p; ++x) {
    // Compute the right-hand side of the elliptic curve equation
    long long squareY = (intPow(x, 3, p) + a * x + b) % p;
    squareY = (squareY + p) % p; // Ensure non-negative result

    // Determine if squareY is a quadratic residue modulo p using the Legendre
    // symbol
    long long legendre = legendreSymbol(squareY, p);

    // If legendre is 1, there are two solutions for y
    if (legendre == 1)
      count += 2;
    // If legendre is 0, there is exactly one solution for y
    else if (legendre == 0)
      count += 1;
  }

  return count; // Return the total number of points on the elliptic curve
}

/**
 * @brief Computes the order of a point on an elliptic curve using the Baby-Step
 * Giant-Step (BSGS) algorithm.
 *
 * @param P The point on the elliptic curve whose order is to be determined.
 * @return The order of the point P.
 */
long long BSGS(Point P) {
  long long p = P.p;
  // Compute Q as P * (p + 1)
  Point Q = P * (p + 1);

  // Check if P is the identity point
  if (P == Point()) {
    return 1;
  }

  // Calculate the step size m
  long long m = static_cast<long long>(std::pow(p, 0.25) + 1);

  // Create a map to store baby steps
  std::unordered_map<Point, long long> babySteps;
  // Initialize temporary point
  Point temp = P;

  // Compute baby steps
  for (long long j = 1; j <= m; ++j) {
    if (temp == Point()) {
      return j;
    }
    babySteps[temp] = j;
    temp = temp + P;
  }

  // Calculate the giant step
  Point giantStep = P * m;
  long long l;
  long long M;

  // Perform giant steps to find collision
  for (long long k = -2 * m; k <= 2 * m; ++k) {
    temp = Q + giantStep * k;
    Point revtemp = temp * (-1);

    // Check for collision in baby steps
    if (babySteps.find(temp) != babySteps.end()) {
      l = babySteps[temp];
      if (p + 1 + m * k - l != 0) {
        M = p + 1 + m * k - l;
        break;
      }
    }
    if (babySteps.find(revtemp) != babySteps.end()) {
      l = babySteps[revtemp];
      if (p + 1 + m * k + l != 0) {
        M = p + 1 + m * k + l;
        break;
      }
    }
  }

exit_loop:
  // Factorize M and check divisibility
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

/**
 * @brief Builds the elliptic curve given by the Weierstrass form y^2 = x^3 + ax
 * + b over the prime field F_p.
 *
 * @param a The coefficient "a" of the elliptic curve equation.
 * @param b The coefficient "b" of the elliptic curve equation.
 * @param p The prime modulus of the elliptic curve.
 *
 * @return A vector of Points representing the points on the elliptic curve.
 */
std::vector<Point> buildCurve(long long a, long long b, long long p) {
  // Initialize the vector of points with the neutral element
  std::vector<Point> points{{0, 0, 0, -1}};

  // Iterate over the range of x coordinates
  for (long long x = -p / 2; x <= p / 2; ++x) {
    // Calculate the square of the y coordinate
    long long squareY = (x * x * x + a * x + b) % p;
    squareY = (squareY + p) % p;

    // Calculate the square root of the y coordinate using the Tonelli-Shanks
    // algorithm
    long long s = sqrtMod(squareY, p);

    // If the square root is not found, skip this x coordinate
    if (s == -1) {
      continue;
    }

    // Calculate the two y roots from the square root
    std::pair<long long, long long> yRoots{s, -s};

    // If the two y roots are the same, add the point with the single y root to
    // the vector of points
    if (yRoots.first == yRoots.second) {
      normalize(yRoots.first, p);
      points.emplace_back(x, yRoots.first, a, p);
      continue;
    }

    // Normalize the y roots
    normalize(yRoots.first, p);
    normalize(yRoots.second, p);

    // Add the two points with the two y roots to the vector of points
    points.emplace_back(x, yRoots.first, a, p);
    points.emplace_back(x, yRoots.second, a, p);
  }

  return points;
}

/// Removes invalid orders from the given set of orders.
/// An order is considered invalid if it is not divisible by the lcmValue.
/// @param lcmValue The least common multiple value used for validation.
/// @param orders A set of orders to be validated and possibly removed.
void removeInvalidOrders(long long lcmValue, std::set<long long> &orders) {
  // Iterate over the set of orders
  for (auto it = orders.begin(); it != orders.end();) {
    // Check if the current order is not divisible by lcmValue
    if (*it % lcmValue != 0) {
      // Erase the current order and update the iterator
      it = orders.erase(it);
    } else {
      // Move to the next order
      ++it;
    }
  }
}

/**
 * @brief Finds the order of the given elliptic curve.
 *
 * @param a The coefficient "a" of the elliptic curve equation.
 * @param b The coefficient "b" of the elliptic curve equation.
 * @param p The prime modulus of the elliptic curve.
 * @return The order of the elliptic curve.
 */
long long findOrder(long long a, long long b, long long p) {
  // Calculate the lower and upper bounds for the possible orders
  long long lowerBound =
      static_cast<long long>(std::ceil(p + 1 - 2 * std::sqrt(p)));
  long long upperBound =
      static_cast<long long>(std::floor(p + 1 + 2 * std::sqrt(p)));

  // Initialize the set of possible orders
  std::set<long long> possibleOrders;

  // Iterate over the range of possible orders and add them to the set
  for (long long N = lowerBound; N <= upperBound; ++N) {
    possibleOrders.insert(N);
  }

  // Initialize the vector of point orders
  std::vector<long long> pointOrders;

  // Seed the random number generator
  std::srand(std::time(0));

  // Iterate until the set of possible orders has only one element
  while (possibleOrders.size() > 1) {
    // Generate a random point on the curve
    long long x = -p / 2 + std::rand() % p;
    long long y = sqrtMod(((x * x * x + a * x + b) % p + p) % p, p);
    if (y > p / 2)
      y -= p;
    if (y == -1)
      continue;
    Point P(x, y, a, p);

    // Calculate the order of the point
    long long pointOrder = BSGS(P);

    // Add the order to the vector of point orders
    pointOrders.push_back(pointOrder);

    // Calculate the least common multiple of the point orders
    long long LCM = findLCM(pointOrders);

    // Remove any orders that are not divisible by the LCM
    removeInvalidOrders(LCM, possibleOrders);
    std::cout << possibleOrders.size() << std::endl;
  }

  // Return the only element in the set of possible orders
  return *possibleOrders.begin();
}

/**
 * @brief Finds all simple order groups of the given elliptic curve.
 *
 * @param a The coefficient "a" of the elliptic curve equation.
 * @param b The coefficient "b" of the elliptic curve equation.
 * @param p The prime modulus of the elliptic curve.
 *
 * @return A vector of vectors of Points, where each inner vector is a simple
 * order group of the elliptic curve.
 */
std::vector<std::vector<Point>> findSimpleOrderGroups(long long a, long long b,
                                                      long long p) {

  std::vector<std::vector<Point>> simpleOrderGroups;
  std::set<Point> gens;

  for (long long x = -p / 2; x < p / 2; ++x) {
    long long squareY = (intPow(x, 3, p) + a * x + b) % p;
    squareY = (squareY + p) % p;
    long long s = sqrtMod(squareY, p);
    if (s == -1) {
      continue;
    }
    normalize(s, p);
    Point point = Point(x, s, a, p);
    bool fl = false;
    long long n = BSGS(point);

    if (isPrimeFermat(n)) {
      std::vector<Point> group;
      Point zero{0, 0, 0, -1};
      Point cur = zero;
      for (long long i = 0; i < n; ++i) {
        group.push_back(cur);
        cur = cur + point;

        if (gens.find(cur) != gens.end()) {
          fl = true;
          break;
        }
      }
      if (!fl) {
        // Add the current point to the set of generators
        gens.insert(point);

        // Add the group to the vector of simple order groups
        simpleOrderGroups.push_back(group);
      }
    }
  }
  // Return the vector of simple order groups
  return simpleOrderGroups;
}
