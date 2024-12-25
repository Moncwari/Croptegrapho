#include "Polynomial.hpp"


// Генерация допустимого образующего полинома для группы
Polynomial polyBuilding(int p, int n) { // p^n
  std::vector<int> coefs(n, 0);
  for (int i{0}; i < n; ++i) {
    int q = std::rand() % p;
    if (i == 0) {
      while (q == 0) {
        q = std::rand() % p;
      }
    }
    coefs[i] = q;
  }
  return Polynomial(coefs, p);
}
// Проверка на неприводимость
bool isIrreducible(Polynomial poly, int p) {
  if (poly.degree() < 4) {
    for (int x = 0; x < p; ++x) {
      int result = 0, power = 1;
      for (int j = 0; j <= poly.degree(); ++j) {
        result = (result + power * poly.getCoef(j)) % p;
        power = (power * x) % p;
      }
      if (result == 0) {
        std::cout << "\nНайден корень " << x << " у многочлена ";
        std::cout << poly;
        return false;
      }
    }
    return true;
  } else {
    for (int i = 0; i < 2; ++i) {
      for (int j = 1; j < (p - 1) * i + 2; ++j) {
        for (int k = 0; k < p; ++k) {
          if (i + j + k != 0) {
            Polynomial square({i, j, k}, p);
            std::cout << square << '\n';
            if (poly % square == Polynomial({0}, p)) {
              std::cout << '\n' << "Поделилось на " << square << '\n';
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

// Вспомогательная функция для итерации
void IncVecOne(std::vector<int> &vecOne, int p, int q) {
  for (auto it = vecOne.rbegin(); it != vecOne.rend(); ++it) {
    int &elem = *it;
    if (elem < q) {
      elem++;
      break;
    }
    elem = 0; // Сбрасываем текущий разряд
  }
}

// Формирование элементов поля
std::vector<Polynomial> fieldBuilding(int n, int p, int q) {
  const int NumIteration = std::pow((q - p + 1), n);
  std::vector<Polynomial> vecRes;
  vecRes.reserve(NumIteration);

  std::vector<int> vecOne(n, 0); // Коэффициенты начинаются с нуля
  for (int i = 0; i < NumIteration; ++i, IncVecOne(vecOne, p, q)) {
    vecRes.emplace_back(Polynomial(vecOne, q));
  }

  return vecRes;
}

// Умножение полиномов в поле
Polynomial fieldMultiply(Polynomial first, Polynomial second,
                         Polynomial formPoly) {
  Polynomial result = first * second;
  result = result % formPoly;
  return result;
}

// Поиск образующих элементов поля
std::vector<Polynomial> findGenerators(const std::vector<Polynomial> &field,
                                       const Polynomial &formPoly, int p) {
  std::vector<Polynomial> generators;
  int fieldOrder = field.size() - 1; // Порядок группы (без нуля)

  for (const auto &elem : field) {
    if (elem == Polynomial({0}, p)) {
      continue; // Пропускаем нулевой элемент
    }
    
    Polynomial cur = elem;
    std::cout << elem << '\n';

    for (int i = 2; i < fieldOrder; ++i)
    {
        cur = fieldMultiply(cur, elem, formPoly);
        cur.trim();
        if (cur == Polynomial({1}, p))
        {
            break;
        }
        
    }

    if (cur == Polynomial({1}, p)) {
      generators.push_back(elem); // Добавляем в список образующих
      std::cout << '\n' << elem << "   !!!\n"; 
    }
  }

  return generators;
}

// Исследование поля
void exploreField(int p, int n, Polynomial formPoly = Polynomial()) {
  if (formPoly.degree() == 0 && formPoly.getCoef(0) == 0) {
    formPoly = polyBuilding(p, n);
    while (!isIrreducible(formPoly, p)) {
      formPoly = polyBuilding(p, n);
    }
    std::cout << "Сгенерирован образующий многочлен: " << formPoly << '\n';
  }
  std::cout << "\n***\n";

  std::vector<Polynomial> field = fieldBuilding(n, 0, p - 1);

  std::cout << "Элементы поля:\n";
  for (Polynomial elem : field) {
    std::cout << elem << '\n';
  }

  std::cout << "\nОбразующие поля:\n";
  std::vector<Polynomial> generators = findGenerators(field, formPoly, p);

  for (Polynomial gen : generators) {
    std::cout << gen << '\n';
  }
}

void printField(const std::vector<Polynomial> &field) {
  for (Polynomial poly : field) {
    std::cout << poly << '\n';
  }
}

int main() {
  exploreField(3, 3);
}