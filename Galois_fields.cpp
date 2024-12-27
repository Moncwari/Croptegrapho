#include "Polynomial.hpp"
#include <bitset>
#include <string>

// Решето Эратосфена для поиска простых делителей
std::vector<int> findPrimeDivisors(int num) {
    std::vector<int> primes;
    std::vector<bool> isPrime(num + 1, true);

    for (int i = 2; i * i <= num; ++i) {
        if (isPrime[i]) {
            for (int j = i * i; j <= num; j += i) {
                isPrime[j] = false;
            }
        }
    }

    for (int i = 2; i <= num; ++i) {
        if (isPrime[i] && num % i == 0) {
            primes.push_back(i);
        }
    }

    return primes;
}

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

// Cложение полиномов в поле
Polynomial fieldSum(Polynomial first, Polynomial second) {
  Polynomial result = first + second;
  result.modulo();
  return result;
}

// Вычитание полиномов поля
Polynomial fieldDif(Polynomial first, Polynomial second) {
  Polynomial result = first - second;
  result.modulo();
  return result;
}

// Умножение полиномов в поле
Polynomial fieldMultiply(Polynomial first, Polynomial second,
                         Polynomial formPoly) {
  Polynomial result = first * second;
  result = result % formPoly;
  return result;
}

// Быстрое возведение в степень
Polynomial fastExponentiation(Polynomial base, int exponent, const Polynomial& modPoly) {
    Polynomial result({1}, base.getDet());
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modPoly;
        }
        base = (base * base) % modPoly;
        exponent /= 2;
    }
    return result;
}

// Проверка, является ли элемент образующим
bool isGeneratorFast(const Polynomial& elem, const Polynomial& modPoly, int p, int fieldOrder) {
    auto primeDivisors = findPrimeDivisors(fieldOrder);
    for (int divisor : primeDivisors) {
        int exponent = fieldOrder / divisor;
        Polynomial power = fastExponentiation(elem, exponent, modPoly);
        if (power == Polynomial({1}, p)) {
            return false; 
        }
    }
    Polynomial finalPower = fastExponentiation(elem, fieldOrder, modPoly);
    return finalPower == Polynomial({1}, p);
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

    for (int i = 1; i < fieldOrder; ++i)
    {
        cur = fieldMultiply(cur, elem, formPoly);
        std::cout << cur;
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

std::vector<Polynomial> findGeneratorsFast(const std::vector<Polynomial>& field, const Polynomial& formPoly, int p) {
    std::vector<Polynomial> generators;
    int fieldOrder = field.size() - 1;

    for (const auto& elem : field) {
        if (elem == Polynomial({0}, p)) continue;
        
        // Проверяем на образующий элемент
        if (isGeneratorFast(elem, formPoly, p, fieldOrder)) {
            generators.push_back(elem);
            std::cout << "Образующий элемент: " << elem << '\n';
        } else {
            std::cout << "Не является образующим: " << elem << '\n';
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

  std::vector<Polynomial> field = fieldBuilding(n, 0, p - 1);

  std::cout << "Элементы поля:\n";
  for (Polynomial elem : field) {
    std::cout << elem << '\n';
  }

  std::cout << "\nОбразующие поля:\n";
  std::vector<Polynomial> generators = findGeneratorsFast(field, formPoly, p);

}

// Разложение по степеням образующего
void degreeDecompose(int p, int n, const Polynomial& formPoly, const Polynomial& genPoly) {
  Polynomial gp = genPoly;
  std::vector<Polynomial> field = fieldBuilding(n, 0, p - 1);
  std::vector<Polynomial> powers;
  powers.push_back(gp);
  for (int i = 0; i < field.size() - 2; ++i) {
    gp = fieldMultiply(gp, gp, formPoly);
    powers.push_back(gp);
  }
  for (const auto& elem : field) {
    if (elem == Polynomial({0}, p)) continue;
    auto pos = std::find(powers.begin(), powers.end(), elem);
    std::cout << elem << " является " << std::distance(powers.begin(), pos) + 1 << " степенью образующего";
  }
}

// Поиск обратного элемента в поле
Polynomial findReverse(const Polynomial& poly, const Polynomial& formPoly, std::vector<Polynomial> field) {
  for (const auto& elem : field) {
    if (elem == Polynomial({0}, formPoly.getDet())) continue;
    if (fieldMultiply(poly, elem, formPoly) == Polynomial({1}, formPoly.getDet())) return elem;
  }
  throw std::invalid_argument("Нет обратного элемента в поле");
}

// Красивый вывод поля
void printField(const std::vector<Polynomial> &field) {
  for (Polynomial poly : field) {
    std::cout << poly << '\n';
  }
}

// Преобразование строки в вектор восьмибитных блоков
std::vector<std::bitset<8>> stringToBinary(const std::string& msg) {
    std::vector<std::bitset<8>> binaryBlocks;
    for (char c : msg) {
        binaryBlocks.emplace_back(static_cast<unsigned long>(c));
    }
    return binaryBlocks;
}

// Преобразование вектора восьмибитных блоков в строку
std::string binaryToString(const std::vector<std::bitset<8>>& binaryBlocks) {
    std::string msg;
    for (const auto& block : binaryBlocks) {
        char c = static_cast<char>(block.to_ulong());
        msg += c;
    }
    return msg;
}

// Преобразование вектора восьмибитных блоков в вектор полиномов 2^8
std::vector<Polynomial> blockTextToPoly(const std::vector<std::bitset<8>>& binaryBlocks) {
  std::vector<Polynomial> polyText;
  for (const auto& block : binaryBlocks) {
    std::vector<int> bits(8);
    for (size_t i = 0; i < 8; ++i) {
        bits[i] = block[i];
    }
    std::reverse(bits.begin(), bits.end());
    Polynomial curPoly = Polynomial(bits, 2);
    polyText.push_back(curPoly);
  }
  return polyText;
}

// Преобразование вектора полиномов в вектор восьмибитных блоков
std::vector<std::bitset<8>> polyBlockTextToBin(const std::vector<Polynomial> polyBlocks) {
  std::vector<std::bitset<8>> blockText;
  for (Polynomial elem : polyBlocks) {
    std::string BSstring;
    for (int i = elem.degree(); i >= 0; --i) {
      BSstring += std::to_string(elem.getCoef(i));
    }
    while (BSstring.size() < 8) {
      BSstring = '0' + BSstring;
    }
    std::bitset<8> BS(BSstring);
    blockText.push_back(BS);
  }
  return blockText;
}

// Зашифрование текста
std::string encryption(std::string opentext, Polynomial a, Polynomial b) {
  Polynomial formPoly = polyBuilding(2, 8);
  std::vector<Polynomial> polyText = blockTextToPoly(stringToBinary(opentext));
  for (auto& elem : polyText) {
    elem = fieldMultiply(elem, a, formPoly);
    elem = fieldSum(elem, b);
  }
  std::string ciphertext = binaryToString(polyBlockTextToBin(polyText));
  return ciphertext;
}

// Расшифрование текста

std::string decryption(std::string ciphertext, Polynomial a, Polynomial b) {
  std::vector<Polynomial> polyText = blockTextToPoly(stringToBinary(ciphertext));
  Polynomial formPoly = polyBuilding(2, 8);
  std::vector<Polynomial> field = fieldBuilding(6, 0, 1);
  Polynomial revA = findReverse(a, formPoly, field);
  for (auto& elem : polyText) {
    elem = fieldDif(elem, b);
    elem = fieldMultiply(elem, revA, formPoly);
  }
  std::string opentext = binaryToString(polyBlockTextToBin(polyText));
  return opentext;
}

int main() {
  std::vector<Polynomial> field = fieldBuilding(6, 0, 1);
  printField(field);
  std::string alphabet;
  for (int i = 0; i <= 255; ++i) {
    char c = static_cast<char>(i);
    alphabet += c;
  }
  std::cout << alphabet << '\n';

  std::vector<std::bitset<8>> bi = stringToBinary(alphabet);
  std::string bis = binaryToString(bi);
  for (auto& elem : bi) {
    std::cout << elem << '\n';
  }
  std::vector<Polynomial> bip = blockTextToPoly(bi);
  for (auto& elem : bip) {
    std::cout << elem << '\n';
  }
  std::vector<std::bitset<8>> bi1 = polyBlockTextToBin(bip);
  for (auto& elem : bi1) {
    std::cout << elem << '\n';
  }
}