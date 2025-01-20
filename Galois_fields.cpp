#include "Polynomial.hpp"
#include <bitset>
#include <sstream>
#include <iterator>
constexpr size_t blockSize = 8;

// Проверка числа на простоту
bool isPrime(int n) {
    for (int i = 2; i*i <= n; i++)
        if (n%i == 0) return false;
    return true;
}

// Развернуть битсет
template <size_t N>
std::bitset<N> reverseBitset(const std::bitset<N>& bits) {
    std::bitset<N> reversed;
    for (size_t i = 0; i < N; ++i) {
        reversed[N - 1 - i] = bits[i];
    }
    return reversed;
}

// Ввод вектора из строки, разделённой пробелами
std::vector<int> getVector() {
    std::vector<int> numbers;
    std::string line;
    while (true) {
        std::cout << "Введите числа через пробел: ";
        std::getline(std::cin, line);
        if (std::cin.fail()) {
            std::cin.clear(); 
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
            std::cout << "Ошибка ввода! Попробуйте снова.\n";
            continue;
        }
        std::istringstream iss(line);
        int num;
        numbers.clear();
        bool error = false;
        while (iss >> num) {
            numbers.push_back(num);
        }
        if (iss.fail() && !iss.eof()) {
            std::cout << "Ошибка: ввод содержит недопустимые символы!\n";
            continue; 
        }
        break;
    }
    return numbers;
}


// Красивый вывод поля
void printField(const std::vector<Polynomial> &field) {
  for (Polynomial poly : field) {
    std::cout << poly << '\n';
  }
}

// Генерация случайного полинома
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
            if (poly % square == Polynomial({0}, p)) {
              return false;
            }
          }
        }
      }
    }
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
    elem = 0;
  }
}

// Формирование элементов поля
std::vector<Polynomial> fieldBuilding(int p, int n) {

  const int NumIteration = std::pow(p, n);
  std::vector<Polynomial> vecRes;
  vecRes.reserve(NumIteration);

  std::vector<int> vecOne(n, 0); // Коэффициенты начинаются с нуля
  for (int i = 0; i < NumIteration; ++i, IncVecOne(vecOne, 0, p - 1)) {
    vecRes.emplace_back(Polynomial(vecOne, p));
  }

  return vecRes;
}

//Получить неприводимый многочлен для построения поля
Polynomial getIrr(int p, int n){
  while (true){
    Polynomial q = polyBuilding(p, n + 1);
    if (isIrreducible(q, p))
      return q;
  }
}

//Вывести все неприводимые полиномы над которыми можно строить поле
std::vector<Polynomial> findIrreducables(int p, int n) {
  std::vector<Polynomial> field = fieldBuilding(p, n + 1);
  std::vector<Polynomial> irreducables;
  for (auto &elem : field) {
    if (isIrreducible(elem, p) && elem.degree() == n) irreducables.push_back(elem);
  }
  return irreducables;
}

// Умножение полиномов в поле
Polynomial fieldMultiply(Polynomial first, Polynomial second,
                         Polynomial formPoly) {
  Polynomial result = first * second;
  result = result % formPoly;
  result.trim();
  return result;
}

// Поиск образующих элементов поля
std::vector<Polynomial> findGenerators(const std::vector<Polynomial> &field,
                                       const Polynomial &formPoly, int p) {
  std::vector<Polynomial> generators = {};
  std::vector<Polynomial> group(field.begin() + 1, field.end());
  int fieldOrder = field.size() - 1; // Порядок группы (без нуля)
  Polynomial zero = Polynomial({0}, p);
  Polynomial singular = Polynomial({1}, p);
  for (const auto &elem : field) {
    if (elem == zero || elem == singular) {
      continue; // Пропускаем нулевой элемент
    }
    bool flag = false;

    Polynomial cur = elem;

    for (int i = 2; i < fieldOrder + 1; ++i) {
      cur = fieldMultiply(cur, elem, formPoly);
      cur.trim();
      if (cur == singular) {
        if (i == fieldOrder) {
          flag = true;
        }
        break;
      }
    }
    if (flag) {
      generators.push_back(elem); // Добавляем в список образующих
    }
  }
  return generators;
}

// Разложение по степеням образующего
void degreeDecompose(int p, int n, const Polynomial &formPoly,
                     const Polynomial &genPoly) {
  Polynomial gp = genPoly;
  std::vector<Polynomial> field = fieldBuilding(p, n);
  std::vector<Polynomial> powers;
  powers.push_back(gp);
  for (int i = 0; i < field.size() - 2; ++i) {
    gp = fieldMultiply(gp, genPoly, formPoly);
    powers.push_back(gp);
  }
  for (const auto &elem : field) {
    if (elem == Polynomial({0}, p))
      continue;
    auto pos = std::find(powers.begin(), powers.end(), elem);
    std::cout << elem << " является " << std::distance(powers.begin(), pos) + 1
              << " степенью образующего\n";
  }
}

// Поиск обратного элемента в поле
Polynomial findReverse(const Polynomial &poly, const Polynomial &formPoly,
                       std::vector<Polynomial> field) {
  Polynomial singular = Polynomial({1}, poly.getDet());
  for (const auto &elem : field) {
    if (fieldMultiply(poly, elem, formPoly) ==
        singular) {
      return elem;
      }
  }
  throw std::invalid_argument("Нет обратного элемента в поле");
}

// Преобразование строки в вектор восьмибитных блоков
std::vector<std::bitset<blockSize>> stringToBinary(const std::string &msg) {
  std::string binaryMsg;
  for (char c : msg) {
    binaryMsg += std::bitset<8>(c).to_string();
  }
  binaryMsg.append((blockSize - binaryMsg.size() % blockSize) % blockSize, '0');
  std::vector<std::bitset<blockSize>> binaryBlocks;
  for (int i = 0; i < binaryMsg.size() ; i += blockSize) {
    binaryBlocks.push_back(std::bitset<blockSize>(binaryMsg.substr(i, blockSize)));
  }
  return binaryBlocks;
}

// Преобразование вектора восьмибитных блоков в строку
std::string binaryToString(const std::vector<std::bitset<blockSize>> &binaryBlocks) {
  std::string binaryMsg;
  for (auto block : binaryBlocks) {
    binaryMsg += block.to_string();
  }
  binaryMsg = binaryMsg.substr(0, binaryMsg.size() / 8 * 8);
  std::string msg;
  for (size_t i = 0; i < binaryMsg.size(); i += 8) {
    msg += static_cast<char>(std::bitset<8>(binaryMsg.substr(i, 8)).to_ulong());
  }
  return msg;
}

// Преобразование вектора восьмибитных блоков в вектор полиномов 2^blockSize
std::vector<Polynomial>
blockTextToPoly(const std::vector<std::bitset<blockSize>> &binaryBlocks) {
  std::vector<Polynomial> polyText;
  for (const auto &block : binaryBlocks) {
    std::vector<int> bits(blockSize);
    for (size_t i = 0; i < blockSize; ++i) {
      bits[i] = block[i];
    }
    std::reverse(bits.begin(), bits.end());
    Polynomial curPoly = Polynomial(bits, 2);
    polyText.push_back(curPoly);
  }
  return polyText;
}

// Преобразование вектора полиномов в вектор восьмибитных блоков
std::vector<std::bitset<blockSize>>
polyBlockTextToBin(const std::vector<Polynomial> polyBlocks) {
  std::vector<std::bitset<blockSize>> blockText;
  for (Polynomial elem : polyBlocks) {
    std::string BSstring;
    for (int i = elem.degree(); i >= 0; --i) {
      BSstring += std::to_string(elem.getCoef(i));
    }
    while (BSstring.size() < blockSize) {
      BSstring = '0' + BSstring;
    }
    std::bitset<blockSize> BS(BSstring);
    blockText.push_back(BS);
  }
  return blockText;
}

// Зашифрование текста
std::string encryption(std::string opentext, Polynomial a, Polynomial b, Polynomial irr) {
  Polynomial formPoly = irr;
  std::vector<Polynomial> polyText = blockTextToPoly(stringToBinary(opentext));
  for (auto &elem : polyText) {
    elem = fieldMultiply(elem, a, formPoly);
    elem = elem + b;
  }
  std::vector<std::bitset<blockSize>> bsvector = polyBlockTextToBin(polyText);

  std::string ciphertext = binaryToString(bsvector);
  return ciphertext;
}

// Расшифрование текста
std::string decryption(std::string ciphertext, Polynomial a, Polynomial b, Polynomial irr) {
  Polynomial formPoly = irr;
  std::vector<std::bitset<blockSize>> bsvector = stringToBinary(ciphertext);
  std::vector<Polynomial> polyText =
      blockTextToPoly(bsvector);

  std::vector<Polynomial> field = fieldBuilding(2, blockSize);
  Polynomial revA = findReverse(a, formPoly, field);

  for (auto &elem : polyText) {
    elem = elem - b;
    elem = fieldMultiply(elem, revA, formPoly);
  }

  std::string opentext = binaryToString(polyBlockTextToBin(polyText));
  return opentext;
}

Polynomial irreducableChoosing(int p, int n) {
  int mode;
  std::cout << "Для операций в поле необходим неприводимый многочлен. Введите 0, если готовы ввести его самостоятельно.\nВведите 1 если нужно сгенерировать его за Вас. Введите 2 если нужно вывести список всех допустимых\nнеприводимых многочленов для данного поля:\n";
  std::cin >> mode;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  if (mode == 0) {
    Polynomial formPoly;
    while (true) {
      std::cout << "Введите коэффициеньты многочлена через пробел, от старшей степени к младшей: \n";
      std::string line;
      std::getline(std::cin, line);
      if (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      std::istringstream iss(line);
      std::vector<int> coefficients;
      int num;
      while (iss >> num) {
        coefficients.push_back(num);
      }
      if (!iss.eof()) {
        std::cout << "Ошибка: ввод содержит недопустимые символы! Попробуйте ещё раз!\n";
        continue;
      }
      Polynomial mayFormPoly = Polynomial(coefficients, p);
      if (mayFormPoly.degree() != n) {
        std::cout << "Ой, со степенью Вашего полинома что-то не так! Попробуйте другой.\n";
        continue;
      }
      if (!isIrreducible(mayFormPoly, p)) {
        std::cout << "Ой, кажется, Ваш полином приводим. Так не пойдёт. Попробуйте ещё раз.\n";
        continue;
      }
      formPoly = mayFormPoly;
      break;
    }
    std::cout << "Замечательно! Ваш неприводимый полином:\n";
    std::cout << formPoly << '\n';
    return formPoly;
  }
  if (mode == 1) {
    Polynomial formPoly = getIrr(p, n);
    std::cout << "Замечательно! Ваш неприводимый полином:\n";
    std::cout << formPoly << '\n';
    return formPoly;
  }
  if (mode == 2) {
    std::vector<Polynomial> irrs = findIrreducables(p, n);
    std::cout << "Допустимые неприводимые полиномы:\n";
    printField(irrs);
    std::cout << "Выберите один из них и введите его коэффициенты от старшей степени к младшей: \n";
    Polynomial formPoly;
    while (true) {
      std::cout << "Введите коэффициенты многочлена через пробел, от старшей степени к младшей: \n";
      std::string line;
      std::getline(std::cin, line);
      if (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      std::istringstream iss(line);
      std::vector<int> coefficients;
      int num;
      while (iss >> num) {
        coefficients.push_back(num);
      }
      if (!iss.eof()) {
        std::cout << "Ошибка: ввод содержит недопустимые символы! Попробуйте ещё раз!\n";
        continue;
      }
      Polynomial mayFormPoly = Polynomial(coefficients, p);
      if (mayFormPoly.degree() != n) {
        std::cout << "Ой, со степенью Вашего полинома что-то не так! Попробуйте другой.\n";
        continue;
      }
      if (!isIrreducible(mayFormPoly, p)) {
        std::cout << "Ой, кажется, Ваш полином приводим. Так не пойдёт. Попробуйте ещё раз.\n";
        continue;
      }
      formPoly = mayFormPoly;
      break;
    }
    std::cout << "Замечательно! Ваш неприводимый полином:\n";
    std::cout << formPoly << '\n';
    return formPoly;
  }
  return Polynomial();
}

int main() {
  while (true) {
    int mode;
    exit_loop:
    std::cout << "Выберите режим. Введите 0 если хотите работать с полем Галуа. Введите 1 если хотите выполнить афинное шифрование или расшифрование: \n";
    std::cin >> mode;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (!mode) {
      int p, n;
      std::cout << "Отлично! Давайте построим поле. Введите p и n. Помните, что p должно быть простым! Ваши p и n через пробел: \n";
      std::cin >> p >> n;
      while (!isPrime(p)) {
        std::cout << "Ошибочка вышла! Ваш p не простой. Попробуйте ввести p и n ещё раз: \n";
        std::cin >> p >> n;
      }
      std::cout << "Принято! Работаем в поле GF(" << p << '^' << n << ')' << '\n';
      std::vector<Polynomial> field = fieldBuilding(p, n);
      Polynomial formPoly = irreducableChoosing(p, n);
      std::cout << "Поле готово. Неприводимый полином:\n" << formPoly << '\n';
      std::cout << "Элементы поля: \n";
      printField(field);
      while (true) {
        std::cout << "Введите 0 если хотите поработать с образующими этого поля. Введите 1, если хотите провести некие арифметические операции. Введите 2 если хотите вернуться назад:\n";
        std::cin >> mode;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (mode == 0) {
          std::cout << "Образующие данного поля: \n";
          std::vector<Polynomial> generators = findGenerators(field, formPoly, p);
          printField(generators);
          Polynomial genPoly;
          while (true) {
            std::cout << "Введите через пробел коэффициенты образующего, по которому хотите разложить данное поле: \n";
            std::vector coefficients = getVector();
            Polynomial mayGenPoly = Polynomial(coefficients, p);
            if (std::find(generators.begin(), generators.end(), mayGenPoly) == generators.end()) {
              std::cout << "Упс! Это не образующий. Попробуйте ещё раз:\n";
              continue;
            }
            genPoly = mayGenPoly;
            break;
          }
          degreeDecompose(p, n, formPoly, genPoly);
          continue;
        }
        if (mode == 1) {
          std::cout << "Введите + если хотите сложить 2 полинома, - если вычесть и * - если умножить:\n";
          std::string mode;
          std::cin >> mode;
          std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          while (!(mode == "+" || mode == "-" || mode == "*")) {
            std::cout << "Недопустимая операция. введите ещё раз:\n";
            std::cin >> mode;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          }
          if (mode == "+") {
            std::cout << "Введите первое слагаемое:\n";
            std::vector<int> coefficients = getVector();
            Polynomial Poly1 = Polynomial(coefficients, p);
            std::cout << "Введите второе слагаемое:\n";
            coefficients = getVector();
            Polynomial Poly2 = Polynomial(coefficients, p);
            std::cout << "Сумма этих полиномов:\n";
            Polynomial sum = Poly1 + Poly2;
            std::cout << sum << '\n';
          }
          if (mode == "-") {
            std::cout << "Введите уменьшаемое:\n";
            std::vector<int> coefficients = getVector();
            Polynomial Poly1 = Polynomial(coefficients, p);
            std::cout << "Введите вычитаемое:\n";
            coefficients = getVector();
            Polynomial Poly2 = Polynomial(coefficients, p);
            std::cout << "Разность этих полиномов:\n";
            Polynomial dif = Poly1 - Poly2;
            std::cout << dif << '\n';
          }
          if (mode == "*") {
            std::cout << "Введите первый множитель:\n";
            std::vector<int> coefficients = getVector();
            Polynomial Poly1 = Polynomial(coefficients, p);
            std::cout << "Введите второй множитель:\n";
            coefficients = getVector();
            Polynomial Poly2 = Polynomial(coefficients, p);
            std::cout << "Произведение этих полиномов:\n";
            Polynomial comp = fieldMultiply(Poly1, Poly2, formPoly);
            std::cout << comp << '\n';
          }
          continue;
        }
        if (mode == 2) {
          goto exit_loop;
        }
      }
      continue;
    }
    if (mode) {
      std::cout << "Введите 0 если хотите осуществить зашифрование, 1 - если расшифрование:\n";
      std::cin >> mode;
      if (!mode) {
        std::cout << "Введите текст для шифрования:\n";
        std::string opentext;
        std::cin >> opentext;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Мы шифруем при помощи блоков длины " << blockSize << ". Это значит, что вы должны подобрать ключи степени " << blockSize - 1 << " и неприводимый полином длины " << blockSize << '\n';
        Polynomial a;
        while (true) {
          std::cout << "Введите коэффициенты первого ключа через пробел:\n";
          std::vector<int> coefficients = getVector();
          a = Polynomial(coefficients, 2);
          std::cout << "Итого Ваш полином: " << a << "\n";
          if (a.degree() != blockSize - 1) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        Polynomial b;
        while (true) {
          std::cout << "Введите коэффициенты второго ключа через пробел:\n";
          std::vector<int> coefficients = getVector();
          b = Polynomial(coefficients, 2);
          std::cout << "Итого Ваш полином: " << b << "\n";
          if (b.degree() != blockSize - 1) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        Polynomial irr;
        while (true) {
          std::cout << "Введите коэффициенты неприводимого полинома поля через пробел:\n";
          std::vector<int> coefficients = getVector();
          irr = Polynomial(coefficients, 2);
          if (irr.degree() != blockSize) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        std::string ciphertext = encryption(opentext, a, b, irr);
        std::cout << "Ваш шифртекст:\n" << ciphertext << '\n';
        continue;
      }
      if (mode) {
        std::cout << "Введите текст для расшифрования:\n";
        std::string ciphertext;
        std::cin >> ciphertext;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Мы шифруем при помощи блоков длины " << blockSize << ". Это значит, что вы должны ввести ключи степени " << blockSize - 1 << " и неприводимый полином длины " << blockSize << '\n';
        Polynomial a;
        while (true) {
          std::cout << "Введите коэффициенты первого ключа через пробел:\n";
          std::vector<int> coefficients = getVector();
          a = Polynomial(coefficients, 2);
          std::cout << "Итого Ваш полином: " << a << "\n";
          if (a.degree() != blockSize - 1) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        Polynomial b;
        while (true) {
          std::cout << "Введите коэффициенты второго ключа через пробел:\n";
          std::vector<int> coefficients = getVector();
          b = Polynomial(coefficients, 2);
          std::cout << "Итого Ваш полином: " << b << "\n";
          if (b.degree() != blockSize - 1) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        Polynomial irr;
        while (true) {
          std::cout << "Введите коэффициенты неприводимого полинома поля через пробел:\n";
          std::vector<int> coefficients = getVector();
          irr = Polynomial(coefficients, 2);
          if (irr.degree() != blockSize) {
            std::cout << "Недопустимый ключ попробуйте ещё раз\n";
          } else break;
        }
        std::string opentext = decryption(opentext, a, b, irr);
        std::cout << "Ваш открытый текст " << opentext << '\n';
        continue;
      }
    }
  }
}
