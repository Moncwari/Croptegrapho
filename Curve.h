#pragma once
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

long long findReverse(long long x, long long p) {
  if (x < 0)
    x = (x % p + p) % p;
  long long a = x, b = p, u = 1, v = 0;
  while (b != 0) {
    long long t = a / b;
    a -= t * b;
    std::swap(a, b);
    u -= t * v;
    std::swap(u, v);
  }
  if (a != 1)
    return -1;
  return (u % p + p) % p;
}

long long powerMod(long long base, long long exp, long long p) {
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
  long long ls = powerMod(a, (p - 1) / 2, p);
  return (ls == p - 1) ? -1 : ls;
}

long long gcd(long long a, long long b) {
  while (b != 0) {
    long long temp = b;
    b = a % b;
    a = temp;
  }
  return a;
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
  if (legendreSymbol(a, p) == -1) {
    return -1;
  }
  if (p % 4 == 3) {
    return powerMod(a, (p + 1) / 4, p);
  }
  return sqrtModTonnelliShanks(a, p);
}
struct Point {
  long long x, y, a, p, order;

  Point() : x(0), y(0), a(0), p(-1) {}

  Point(long long _x, long long _y, long long _a, long long _p)
      : x(_x), y(_y), a(_a), p(_p) {}

  static void normalize(long long &x, long long p) {
    x %= p;
    if (x > (p - 1) / 2) {
      x -= p;
    }
    if (x < -(p - 1) / 2) {
      x += p;
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
      long long rev_y = findReverse(2 * y, p);
      resX = (powerMod(((3 * x * x + a) % p) * rev_y, 2, p) - 2 * x);
      normalize(resX, p);
      resY = (((((3 * x * x + a) % p) * rev_y) % p) * (x - resX) - y);
      normalize(resY, p);
    } else {
      long long rev_y = findReverse(addend.x - x, p);
      normalize(rev_y, p);
      resX = powerMod((addend.y - y) * rev_y, 2, p) - x - addend.x;
      normalize(resX, p);
      resY = ((addend.y - y) * rev_y) % p * (x - resX) - y;
      normalize(resY, p);
    }
    return Point(resX, resY, a, p);
  }

  Point operator*(long long scalar) const {
    Point base = *this;
    if (scalar == 0 || p == -1)
      return Point();
    if (scalar == -1)
      return Point(x, -y, a, p);
    if (scalar < 0)
      return (base * (-scalar)) * (-1);

    Point result{x, y, a, p};
    scalar--;

    while (scalar) {
      if (scalar & 1) {
        result = result + base;
      }
      base = base + base;
      scalar >>= 1;
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

long long nativeFindOrder(long long a, long long b, long long p) {
  long long count = 1; // Start with the neutral element (point at infinity)

  // Iterate over all possible x values in the field F_p
  for (long long x = 0; x < p; ++x) {
    // Compute the right-hand side of the elliptic curve equation
    long long squareY = (powerMod(x, 3, p) + a * x + b) % p;
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

long long BSGS(Point P) {
  long long p = P.p;
  long long m = static_cast<long long>(std::pow(p, 0.25) + 1);

  Point Q = P * (p + 1);
  if (P == Point()) {
    return 1;
  }
  std::unordered_map<Point, long long> babySteps;
  Point temp = P;

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

long long findOrder(long long a, long long b, long long p) {
  long long lowerBound =
      static_cast<long long>(std::ceil(p + 1 - 2 * std::sqrt(p)));
  long long upperBound =
      static_cast<long long>(std::floor(p + 1 + 2 * std::sqrt(p)));
  std::vector<long long> possibleOrders;

  for (long long N = lowerBound; N <= upperBound; ++N) {
    possibleOrders.push_back(N);
  }
  long long ans = 0;
  long long LCM = 1;
  std::vector<long long> pointOrders;
  for (long long x = -(p - 1) / 2; x < (p + 1) / 2; ++x) {
    long long y = sqrtMod((powerMod(x, 3, p) + a * x + b), p);
    if (y == -1)
      continue;
    normalize(y, p);
    Point P(x, y, a, p);
    long long s = BSGS(P);
    LCM = lcm(LCM, s);
    long long cou = 0;
    for (long long i = lowerBound; i <= upperBound; ++i) {
      if (i % LCM == 0) {
        cou += 1;
        ans = i;
      }
    }
    if (cou == 1) {
      return ans;
    }
  }
  return LCM;
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
    long long squareY = (powerMod(x, 3, p) + a * x + b) % p;
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

Point randPoint(long long a, long long b, long long p) {
  while (true) {
    long long x = -p / 2 + std::rand() % p;
    long long y = sqrtMod((powerMod(x, 3, p) + a * x + b), p);
    normalize(y, p);
    if (y == -1)
      continue;
    Point P(x, y, a, p);
    return P;
  }
}

std::vector<std::vector<Point>> findSimpleOrderGroups1(long long a, long long b,
                                                       long long p) {
  long long order = findOrder(a, b, p);
  Point P = randPoint(a, b, p);
  while (BSGS(P) != order) {
    P = randPoint(a, b, p);
  }
  std::vector<long long> factors = factorize(order);
  std::vector<std::vector<Point>> groups;
  std::unordered_map<int, int> freq;
  for (long long f : factors)
    freq[f]++;
  std::vector<int> ufactors;
  for (int x : factors) {
    if (freq[x] == 1)
      ufactors.push_back(x);
  }
  for (int x : ufactors) {
    Point gen = P * (order / x);
    Point p = gen;
    std::vector<Point> subgroup;
    subgroup.push_back(Point());
    do {
      subgroup.push_back(p);
      p = p + gen;
    } while (!(p == Point()));
    groups.push_back(subgroup);
  }
  return groups;
}