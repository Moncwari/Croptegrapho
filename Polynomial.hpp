#pragma once

#include <algorithm>
#include <climits>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <set>
class Polynomial {

private:
  std::vector<int> coefficients; // От младшей степени к старшей внутри
                                 // реализации, от старшей к младшей при вводе
  int galuaDet;

public:
  int getDet() const { return galuaDet; }

  int getCoef(const int degree) const { return coefficients[degree]; }

  // Функция удаления ведущих нулей в массиве коэффициентов
  void trim() {
    while (coefficients.size() > 1 && coefficients.back() == 0) {
      coefficients.erase(coefficients.end());
    }
  }
  // Приведение многочлена по модулю
  void modulo() {
    for (size_t i = 0; i < coefficients.size(); ++i) {
      coefficients[i] = (coefficients[i] % galuaDet + galuaDet) % galuaDet;
    }
  }

  // По умолчанию инициализирую многочлен как нулевой
  Polynomial() {
    coefficients = std::vector<int>({0});
    galuaDet = INT_MAX;
  }

  // Стандартный конструктор от вектора коэффициентов
  Polynomial(const std::vector<int> &coeffs, int det = INT_MAX)
      : galuaDet{det},
        coefficients{coeffs.rbegin(), std::reverse_iterator{std::find_if(
                                          coeffs.begin(), coeffs.end() - 1,
                                          [](int a) { return a; })}} {}
  // Получить степень многочлена
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

  // Перегрузка оператора сложения
  Polynomial operator+(const Polynomial &addend) const {
    size_t sumDegree =
        std::max(coefficients.size(), addend.coefficients.size());
    std::vector sum(sumDegree, 0);

    for (size_t i = 0; i < sumDegree; ++i) {
      int a = (i < coefficients.size()) ? coefficients[i]
                                        : 0; // (Условие) ? Если да : Если нет
      int b = (i < addend.coefficients.size()) ? addend.coefficients[i] : 0;
      sum[i] = a + b;
    }

    return Polynomial(sum);
  }

  // Перегрузка оператора вычитания
  Polynomial operator-(const Polynomial &subtrahend) const {
    size_t diffDegree =
        std::max(coefficients.size(), subtrahend.coefficients.size());
    std::vector<int> difference(diffDegree, 0);

    for (size_t i = 0; i < diffDegree; ++i) {
      int a = (i < coefficients.size()) ? coefficients[i]
                                        : 0; // Коэффициент из уменьшаемого
      int b = (i < subtrahend.coefficients.size())
                  ? subtrahend.coefficients[i]
                  : 0; // Коэффициент из вычитаемого
      difference[i] = a - b;
    }

    return Polynomial(difference);
  }

  // Перегрузка оператора умножения
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

    return Polynomial(comp);
  }

  // Перегрузка оператора взятия остатка от деления
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

  // Перегрузка оператора равенства
  bool operator==(const Polynomial &other) const {
    return (coefficients == other.coefficients && galuaDet == other.galuaDet);
  }

  // Перегрузка вывода полинома
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
        // Добавляем знак для всех членов, кроме первого
        if (!isFirst) {
          os << (coeff > 0 ? " + " : " - ");
        } else if (coeff < 0) {
          os << "-";
        }

        // Вывод коэффициента
        if (std::abs(coeff) != 1 || currentDegree == 0) {
          os << std::abs(coeff);
        }

        // Вывод переменной и степени
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

    // Если все коэффициенты равны 0
    if (isFirst) {
      os << "0";
    }

    return os;
  }
};