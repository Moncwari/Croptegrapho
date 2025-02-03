from itertools import product
import math
import sympy
import random
from numpy import sign

class EllCurve:
    def __init__(self, p, a, b):
        self.p = p
        self.a = a
        self.b = b

    def square_root(self, k): # Корень из k по модулю p
        self.k = k
        if sympy.legendre_symbol(self.k, self.p) == -1:
            return 0.5
        else:
            if sympy.legendre_symbol(-1, self.p) == -1:
                b = -1
            else:
                b = 2
                while sympy.legendre_symbol(b, self.p) != -1:
                    b += 1
            s = 0
            p_1 = self.p-1
            while p_1 % 2 == 0:
                p_1 = p_1 // 2
                s += 1
            t = (self.p-1) // 2**s

            reversed_k = rev_a(self.p, self.k)
            C_0 = pow(b, t, self.p)
            r = pow(self.k, (t+1)//2, self.p)

            for i in range(1, s):
                d_i = pow((r**2*reversed_k), 2**(s-i-1), self.p)
                if d_i % self.p == self.p-1:
                    r = (r*C_0) % self.p
                C_0 = pow(C_0, 2, self.p)
            return to_sym_norm_F(r, self.p, 0), -to_sym_norm_F(r, self.p, 0)

    def summary(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        #print(self.point1, self.point2)
        x1, y1 = self.point1[0], self.point1[1]
        x2, y2 = self.point2[0], self.point2[1]
        if self.point1 == ("a", "a"): # 0 + P = P
            return self.point2
        elif self.point2 == ("a", "a"): # P + 0 = P
            return self.point1
        elif y1 == 0 and y2 == 0: # если P = (a, 0), то 2P = 0
            return ("a", "a")
        else:
            if x1 != x2: # P != Q, x1 != x2
                rev_div = rev_a(self.p, to_sym_norm_F((x2-x1), self.p, 1))
                x3 = ((y2-y1)*rev_div)**2 - x1 - x2
                x3 = to_sym_norm_F(x3, self.p, 0)
                y3 = (y2-y1)*rev_div*(x1-x3) - y1
                y3 = to_sym_norm_F(y3, self.p, 0)
                return (x3, y3)
            else:
                if y1 == y2: # P = Q, x1 = x2, y1 = y2
                    rev_div = rev_a(self.p, (2*y1)%self.p)
                    x3 = ((3*x1**2+self.a)*rev_div)**2 - 2*x1
                    x3 = to_sym_norm_F(x3, self.p, 0)
                    y3 = ((3*x1**2+self.a)*rev_div)*(x1-x3) - y1
                    y3 = to_sym_norm_F(y3, self.p, 0)
                    return (x3, y3)
                else: # P =- Q, x1 = x2, y1 = -y2
                    return ("a", "a")

    def minus_point(self, point):
        self.point = point
        if self.point == ("a", "a"):
            return ("a", "a")
        else:
            return (self.point[0], -self.point[1])

    def point_degree(self, point, degree):
        self.point = point
        self.degree = degree
        flag = False
        if self.degree < 0:
            flag = True
        self.degree = abs(degree)
        if self.degree == 0:
            return ("a", "a")

        degrees = []
        while self.degree > 0:
            i = 1
            while i*2 <= self.degree:
                i *= 2
            else:
                self.degree -= i
                degrees.append(i)

        deg = degrees[0]
        new_point = self.point
        counter = 1
        while counter < deg:
            new_point = self.summary(new_point, new_point)
            counter *= 2
        ans = new_point
        if len(degrees) == 1:
            if flag:
                return self.minus_point(ans)
            else:
                return ans
        else:
            for i in range (1, len(degrees)):
                new_point = self.point
                counter = 1
                while counter < degrees[i]:
                    new_point = self.summary(new_point, new_point)
                    counter *= 2
                ans = self.summary(new_point, ans)
            if flag:
                return self.minus_point(ans)
            else:
                return ans

    def point_creator(self, x):
        self.x = x
        sqr_y = (x**3 + self.a*x + self.b) % self.p
        if sqr_y == 0:
            return ((x, 0), (x, 0))
        else:
            y = self.square_root(sqr_y)
            if y != 0.5:
                return (x, y[0]), (x, y[1])
            else:
                return -1

    def baby_giant(self, point):
        self.P = point
        Q = self.point_degree(self.P, (self.p+1))
        #print(Q)
        m = int((self.p)**(1/4)) + 1
        jPs = [self.point_degree(self.P, j) for j in range(0, m+1)]

        for k in range (-m, m+1):
            Q_ = self.summary(Q, self.point_degree(self.P, 2*m*k))
            if Q_ in jPs:
                l = -jPs.index(Q_)
                M = self.p + 1 + 2*m*k + l
                if M == 0:
                    continue
                else:
                    break
            else:
                if self.minus_point(Q_) in jPs:
                    l = jPs.index(self.minus_point(Q_))
                    M = self.p + 1 + 2*m*k + l
                    if M == 0:
                        continue
                    else:
                        break

        M = self.p + 1 + 2*m*k + l
        fact_M = factor(M)
        for p_i in fact_M:
            new_M = M//p_i
            if self.point_degree(self.P, new_M) == ("a", "a"):
                M = new_M
        return M

    def order(self):
        LCM = 1
        flag = False
        for x in range (-(self.p-1)//2, (self.p+1)//2):
            point = self.point_creator(x)
            if point != -1:
                LCM = math.lcm(LCM, self.baby_giant(point[0]))
                counter = 0
                for i in range (int(self.p+1-2*self.p**(1/2))+1, int(self.p+1+2*self.p**(1/2))+1):
                    if i % LCM == 0:
                        N = i
                        counter += 1
                if counter == 1:
                    return N

    def find_subgropus(self):
        N = self.order()
        p_0 = ("a", "a")
        fact_N = factor(N)
        list_fact_N = list(factor(N))
        point_in_subgroups = []
        for _ in range (len(fact_N)):
            point_in_subgroups.append(set([p_0]))
        subgropus_dict = dict(zip(list_fact_N, point_in_subgroups))

        for x in range (-(self.p-1)//2, (self.p+1)//2):
            point = self.point_creator(x)
            #print(point)
            if point != -1:
                point_order = self.baby_giant(point[0])
                if (point_order in fact_N) and (point[0] not in subgropus_dict[point_order]) and (point[1] not in subgropus_dict[point_order]): 
                    subgropus_dict[point_order].add(point[0])
                    subgropus_dict[point_order].add(point[1])
                    print("Подгруппа порядка", point_order)
                    for i in range (point_order):
                        point_degr = self.point_degree(point[0], i)
                        print(point_degr)
                        subgropus_dict[point_order].add(point_degr)
                    print()
        return subgropus_dict

def rev_a(n, a):
    r, y1, y2 = 1, 1, 0
    n_old = n
    while r:
        q = n//a
        r = n - q*a
        y = y2 - q*y1
        n = a
        a = r
        y2 = y1
        y1 = y
    return (y2 + n_old) % n_old

def to_sym_norm_F(a, p, key):
    if key == 0: # В симметричное поле
        a = a%p
        if a > (p-1)//2:
            a -= p
    if key == 1: # В обычное поле
        a = ((a%p) + p) % p
    return a

#print(to_sym_norm_F(square_root(10, 13)[0], 13, 0), to_sym_norm_F(square_root(10, 13)[1], 13, 0))

def factor(a):
    ans = set()
    if a % 2 == 0:
        ans.add(2)
        while a % 2 == 0:
            a //= 2
    d = 3
    while a > 1:
        if a % d == 0:
            ans.add(d)
            while a % d == 0:
                a //= d
        else:
            d += 2
    return ans

def Ferma(n):
    flag = True
    for t in range (1, 10):
        a = random.randint(2, n-1)
        r = pow(a, n-1, n)
        if r != 1:
            flag = False
    if flag:
        return n
    else:
        return 0
#print(Ferma(63973))

#print(rev_a(5, -2))
#print(factor(13))

#for n in range (2**45, 2**45+1000):
    #if Ferma(n):
        #print(n)
curve = EllCurve(997, -2, 1)
print("Порядок:", curve.order())
print(curve.find_subgropus())
        #break

#curve = EllCurve(13, 2, 2)
#for x in range(-6, 7):
    #point = curve.point_creator(x)
    #if point != -1:
        #print(point[0])
        #print(curve.baby_giant(point[0]))

print("""Введите 1, если хотите построить эллиптическую кривую,
      отобразить группу её точек и вычислить порядок таковой""")
print("Введите 2, если хотите вычислить точку заданной кратности")
print("""Введите 3, если хотите найти все подгруппы группы точек
      эллиптической кривой некоторого порядка""")
a = int(input())
if a == 1:
    print("""Введите через пробел: p - порядок простого
          конечного поля, a и b - коэффициенты в уравнении
          эллиптической кривой y^2 = x^3 + ax + b""")
    string = input().split()
    p, a, b = int(string[0]), int(string[1]), int(string[2])
    if (-4*a**3 - 27*b**2) % p == 0:
        print("""Введённые значения не соответствуют условию
              гладкости эллиптической кривой""")
    else:
        curve = EllCurve(p, a, b)
        a = curve.point_creator()
        print(a)
        #b = curve.point_degree([-5, -6], 3)
        #print(b)
        #for point in a:
            #print(point, curve.baby_giant(point))
            #print()
