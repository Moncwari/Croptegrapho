#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <initializer_list>
#include <utility>
#include <algorithm>
#include "Polynomial.hpp"



// Генерация допустимого образующего полинома группы
Polynomial polyBuilding(int p, int n)
{ // p^n
    std::vector<int> coefs(n, 0);
    for (int i{0}; i < n; ++i)
    {
        coefs[i] = -p + (rand() % (2 * p + 1));
    }
    return Polynomial(coefs, p);
}

// Проверка на неприводимость
bool isIrreducible(const Polynomial poly, int p)
{
    if (poly.degree() < 4)
    {
        for (int x = 0; x < p; ++x)
        {
            int result = 0, power = 1;
            for (int j = 0; j <= poly.degree(); ++j)
            {
                result = (result + power * poly.getCoef(j)) % p;
                power = (power * x) % p;
            }
            if (result == 0)
            {
                std::cout << '\n'
                          << "Найден корень " << x << '\n';
                return false;
            }
        }
        return true;
    }
    else
    {
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 1; j < (p - 1) * i + 2; ++j)
            {
                for (int k = 0; k < p; ++k)
                {
                    if (i + j + k != 0)
                    {
                        Polynomial square({i, j, k}, p);
                        std::cout << square << '\n';
                        if (poly % square == Polynomial({0}, p))
                        {
                            std::cout << '\n'
                                      << "Поделилось на " << square << '\n';
                            return false;
                        }
                    }
                }
            }
        }
        std::cout << "It's irreducible.\n";
        return true;
    }
}

void exploreGroup(int p, int n, Polynomial formPoly = Polynomial())
{
    if (formPoly.degree() == 0 && formPoly.getCoef(0) == 0)
    {
        formPoly = polyBuilding(p, n);
        while (!isIrreducible(formPoly, p))
        {
            formPoly = polyBuilding(p, n);
        }
        std::cout << "Для Вашего поля был сгенерирован образующий многочлен: " << formPoly;
    }
}

int main()
{
    Polynomial a = Polynomial({11, 4, 6, 5, 2}, 11);
    Polynomial b = Polynomial();
    Polynomial c = a % b;
    std::cout << c;
}