from itertools import product
import math
import sympy

class EllCurve:
    def __init__(self, p, a, b):
        self.p = p
        self.a = a
        self.b = b

    def summary(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        if self.point1 == [0, 0]:
            return self.point1
        elif self.point2 == [0, 0]:
            return self.point1
        else:
            if self.point1[0] != self.point2[0]: # P!=Q, x1!=x2
                rev_div = rev_a(self.p, (point2[0]-point1[0])%self.p)
                x3 = ((point2[1]-point1[1])*rev_div)**2 - point1[0] - point2[0]
                y3 = ((point2[1]-point1[1])*rev_div) * (point1[0] - x3%self.p) - point1[1]
                x3, y3 = x3%self.p, y3%self.p
                if x3 < -(self.p-1)//2 or (self.p-1)//2 < x3:
                    x3 -= self.p*(x3//abs(x3))
                if y3 < -(self.p-1)//2 or (self.p-1)//2 < y3:
                    y3 -= self.p*(y3//abs(y3))
                return [x3, y3]
            else:
                if self.point1[0] == self.point2[0] and self.point1[1] == self.point2[1]: # P=Q
                    rev_div = rev_a(self.p, (2*point1[1])%self.p)
                    x3 = ((3*point1[0]**2+self.a)*rev_div)**2 - 2*point1[0]
                    y3 = ((3*point1[0]**2+self.a)*rev_div) * (point1[0] - x3%self.p) - point1[1]
                    x3, y3 = x3%self.p, y3%self.p
                    if x3 < -(self.p-1)//2 or (self.p-1)//2 < x3:
                        x3 -= self.p*(x3//abs(x3))
                    if y3 < -(self.p-1)//2 or (self.p-1)//2 < y3:
                        y3 -= self.p*(y3//abs(y3))
                    return [x3, y3]
                else: # P=-Q, x1=x2, y1=-y2
                    return [0, 0]

    def point_degree(self, point, degree):
        self.point = point
        self.degree = degree
        flag = False
        if self.degree < 0:
            flag = True
        self.degree = abs(degree)
        if self.degree == 0:
            return [0, 0]
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
                return [ans[0], -ans[1]]
            else:
                return ans
        else:
            for i in range (1, len(degrees)):
                counter = 1
                new_point = self.point
                while counter < degrees[i]:
                    new_point = self.summary(new_point, new_point)
                    counter *= 2
                ans = self.summary(new_point, ans)
            if flag:
                return [ans[0], -ans[1]]
            else:
                return ans

    def point_creator(self):
        squares = set()
        half_F_squares = []
        points = []
        for i in range ((self.p+1)//2):
            sqr = i**2 % self.p
            squares.add(sqr)
            half_F_squares.append(sqr)
        print(squares)
        for x in range (-(self.p-1)//2, (self.p+1)//2):
            sqr_y = (x**3 + self.a*x + self.b) % p
            if sqr_y in squares:
                index = half_F_squares.index(sqr_y)
                if index == 0:
                    points.append([x, index])
                else:
                    points.append([x, index])
                    points.append([x, -(index)])
        return points

    def baby_giant(self, point):
        self.point = point
        Q = self.point_degree(self.point, (self.p+1))
        m = int(self.p**(1/4)) + 1
        print(m)
        jPs = []
        for j in range (0, m+1):
            jPs.append(self.point_degree(self.point, j))
        print(jPs)
        for k in range (-m, m+1):
            P_ = self.point_degree(self.point, 2*m*k)
            Q_ = self.summary(Q, P_)
            print(k, Q_)
            flag = False
            for jP in jPs:
                if jP[0] == Q_[0] and abs(jP[1]) == abs(Q_[1]):
                    if jP == [0, 0]:
                        l = -1
                        print("l =", l)
                        flag = True
                        break
                    else:
                        l = -Q_[1]//jP[1]
                        print("l =", l)
                        flag = True
                        break
            if flag:
                print("!!!")
                break
            
        flag1 = False
        M = self.p + 1 + 2*m*k + l
        print(M)
        while not flag1:
            flag2 = False
            factor_M = list(factor(M))
            for p in factor_M:
                P_ = self.point_degree(self.point, M//p)
                if P_ == [0, 0]:
                    M = M//p
                    flag2 = True
            if not flag2:
                flag1 = True
        return M

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

def factor(a):
    ans = set()
    d = 2
    while a > 1:
        if sympy.isprime(d) and a%d==0:
            ans.add(d)
            a //= d
        else:
            d += 1
    return ans


#print(rev_a(267, 13))
#print(factor(13))
#a = curve.curve_creator()
#print(a)

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
        for point in a:
            print(point, curve.baby_giant(point))
            print()
