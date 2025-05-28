#ifndef BIGINT_HPP
#define BIGINT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <bitset>

class BigInt {
    public:
        std::vector<uint64_t> digits;
        static const uint64_t BASE = 1'000'000'000'000'000'000ULL;
        bool is_negative;

        BigInt() : is_negative(false) {}

        BigInt(int64_t num) {
            is_negative = (num < 0);
            if (is_negative) num = -num;
            while (num) {
                digits.push_back(num % BASE);
                num /= BASE;
            }
            if (digits.empty()) digits.push_back(0);
        }

        BigInt(const std::vector<uint64_t>& d, bool neg = false) : digits(d), is_negative(neg) {
            remove_leading_zeros();
        }

        BigInt(const std::string &s) {
            is_negative = (s[0] == '-');
            for (int i = s.size(); i > is_negative; i -= 18) {
                int len = (i >= 18 + is_negative) ? 18 : i - is_negative;
                digits.push_back(std::stoull(s.substr(i - len, len)));
            }
            remove_leading_zeros();
        }

        BigInt& operator=(const BigInt& other) {
            if (this != &other) {
                digits = other.digits;
                is_negative = other.is_negative; 
            }
            return *this;
        }
        
        bool operator==(const BigInt &other) const {
            return is_negative == other.is_negative && digits == other.digits;
        }

        BigInt operator+(const BigInt &other) const {
            if (is_negative == other.is_negative) {
                BigInt res = add_abs(*this, other);
                res.is_negative = is_negative;
                return res;
            }
            return (*this < other) ? sub_abs(other, *this).negate() : sub_abs(*this, other);
        }
        
        BigInt operator-(const BigInt &other) const {
            if (is_negative != other.is_negative) {
                BigInt res = add_abs(*this, other);
                res.is_negative = is_negative;
                return res;
            }
            return (*this < other) ? sub_abs(other, *this).negate() : sub_abs(*this, other);
        }

        BigInt operator-() const {
            BigInt result = *this;
            if (result != BigInt(0)) {
                result.is_negative = !result.is_negative;
            }
            return result;
        }
    
        friend std::ostream &operator<<(std::ostream &out, const BigInt &num) {
            if (num.is_negative) out << "-";
            out << num.digits.back();
            for (int i = num.digits.size() - 2; i >= 0; --i) {
                out.width(18);
                out.fill('0');
                out << num.digits[i];
            }
            return out;
        }
    
        bool operator<(const BigInt &other) const {
            if (is_negative != other.is_negative) return is_negative;
            if (digits.size() != other.digits.size())
                return (digits.size() < other.digits.size()) ^ is_negative;
            for (int i = digits.size() - 1; i >= 0; --i)
                if (digits[i] != other.digits[i])
                    return (digits[i] < other.digits[i]) ^ is_negative;
            return false;
        }

        bool operator!=(const BigInt &other) const {
            return !(*this == other);
        }
        
        bool operator<=(const BigInt &other) const {
            return (*this < other) || (*this == other);
        }
        
        bool operator>(const BigInt &other) const {
            return !(*this <= other);
        }
        
        bool operator>=(const BigInt &other) const {
            return !(*this < other);
        }
        
        BigInt operator*(const BigInt &other) const {
            BigInt res = karatsuba_mult(*this, other);
            res.is_negative = (this->is_negative != other.is_negative);
            if (res.digits.size() == 1 && res.digits[0] == 0)
                res.is_negative = false;
            return res;
        }
        
        BigInt operator*(const int &oth) const {
            BigInt other = BigInt(oth);
            BigInt res = karatsuba_mult(*this, other);
            res.is_negative = (this->is_negative != other.is_negative);
            if (res.digits.size() == 1 && res.digits[0] == 0)
                res.is_negative = false;
            return res;
        }

        BigInt operator/(const BigInt &other) const {
            if (other == BigInt(0))
                throw std::runtime_error("Division by zero");
        
            BigInt dividend = *this, divisor = other;
            dividend.is_negative = divisor.is_negative = false;
        
            BigInt quotient, current;
            quotient.digits.resize(dividend.digits.size());
        
            for (int i = dividend.digits.size() - 1; i >= 0; --i) {
                current = current.shift(1);
                current.digits[0] = dividend.digits[i];
                current.remove_leading_zeros();
        
                uint64_t x = 0, l = 0, r = BASE;
                while (l <= r) {
                    uint64_t m = (l + r) / 2;
                    if (current < mul_abs(divisor, m))
                        r = m - 1;
                    else {
                        x = m;
                        l = m + 1;
                    }
                }
                
                quotient.digits[i] = x;
                current = sub_abs(current, mul_abs(divisor, x));
            }
        
            quotient.is_negative = (is_negative != other.is_negative);
            quotient.remove_leading_zeros();
            return quotient;
        }
        
        BigInt operator%(const BigInt &other) const {
            BigInt res = *this - (*this / other) * other;
            if (res.is_negative) {
                res = res + other;
            }
            return res;
        }

        BigInt operator^(const BigInt &exp) const {
            if(exp.is_negative)
                throw std::runtime_error("Negative exponent not supported");
        
            BigInt result(1);
            BigInt base = *this;
            BigInt exponent = exp;
            BigInt zero(0);
            BigInt two(2);
        
            while (!(exponent == zero)) {
                if (!(exponent % two == BigInt(0))) {
                    result = result * base;
                }
                base = base * base;
                exponent = exponent / two;
            }
            
            return result;
        }

        BigInt operator&(const BigInt& other) const {
            BigInt result;
            size_t max_size = std::max(digits.size(), other.digits.size());
            result.digits.resize(max_size);
        
            for (size_t i = 0; i < max_size; ++i) {
                uint64_t a = (i < digits.size()) ? digits[i] : 0;
                uint64_t b = (i < other.digits.size()) ? other.digits[i] : 0;
                result.digits[i] = a & b;
            }
        
            result.remove_leading_zeros();
            return result;
        }

        BigInt operator<<(int shift) const {
            if (shift < 0) throw std::invalid_argument("Negative shift");
            std::bitset<1024> bits = bigIntToBits(*this);
            bits <<= shift;
            return bitsToBigInt(bits);
        }
    
        BigInt operator>>(int shift) const {
            if (shift < 0) throw std::invalid_argument("Negative shift");
            std::bitset<1024> bits = bigIntToBits(*this);
            bits >>= shift;
            return bitsToBigInt(bits);
        }
    
        static std::bitset<1024> bigIntToBits(const BigInt& num) {
            std::bitset<1024> bits;
            BigInt n = num;
            int i = 0;
            while (n > BigInt(0)) {
                bits[i++] = ((n % BigInt(2)) == BigInt(1));
                n = n / BigInt(2);
            }
            return bits;
        }
    
        static BigInt bitsToBigInt(const std::bitset<1024>& bits) {
            BigInt result(0);
            for (int i = bits.size() - 1; i >= 0; --i) {
                result = result * BigInt(2) + BigInt(static_cast<int>(bits[i]));
            }
            return result;
        }

        std::vector<uint8_t> toBytes() const {
            std::vector<uint64_t> digits = this->digits;
            std::vector<uint8_t> output;
            output.reserve(digits.size() * 8);
            for (uint64_t value : digits) {
                for (int i = 0; i < 8; ++i) {
                    output.push_back(static_cast<uint8_t>(value >> (i * 8)) & 0xFF);
                }
            }
                
            return output;
        }

        BigInt BytesToBigInt(const std::vector<uint8_t>& input) {
            std::vector<uint64_t> output;
            const size_t num_blocks = (input.size() + 7) / 8;
            output.reserve(num_blocks);

            for (size_t i = 0; i < input.size(); i += 8) {
                uint64_t value = 0;
                
                for (int j = 0; j < 8; ++j) {
                    const size_t byte_index = i + j;
                    const uint8_t byte = (byte_index < input.size()) ? input[byte_index] : 0;
                    value |= static_cast<uint64_t>(byte) << (j * 8);
                }
                
                output.push_back(value);
            }

            return BigInt(output);
        }

        BigInt sqrt() const {
            if (is_negative)
                throw std::runtime_error("Square root of negative number");
                
            if (*this == 0 || *this == 1)
                return *this;
                
            BigInt low = 0;
            BigInt high = *this;
            BigInt res = 0;
            BigInt mid, square;
            BigInt two = 2;
            
            while (low <= high) {
                mid = (low + high) / two;
                square = mid * mid;
                
                if (square == *this)
                    return mid;
                    
                if (square < *this) {
                    low = mid + 1;
                    res = mid;
                } else {
                    high = mid - 1;
                }
            }
            
            return res;
        }

        BigInt gcd(BigInt a, BigInt b) {
            while (b != 0) {
                BigInt tmp = b;
                b = a % b;
                a = tmp;
            }
            return a;
        }

        BigInt bitLength() const { 
            BigInt res = 0;
            BigInt value = *this;
            while (value != 0) {
                value = value / 2;
                res = res + 1;
            }
            return res;
        }
                
        BigInt modExp(const BigInt &exponent, const BigInt &mod) const {
            BigInt base = *this;
            BigInt result(1);
            BigInt b = base % mod;
            BigInt e = exponent;
            BigInt zero(0);
            BigInt two(2);
            while (!(e == zero)) {
                if (!(e % two == BigInt(0)))
                    result = (result * b) % mod;
                b = (b * b) % mod;
                e = e / two;
            }
            return result;
        }


        BigInt findReverse(BigInt p) {
            BigInt x = *this;
            if (x < BigInt(0))
                x = (x % p + p) % p;
            BigInt a = x, b = p, u = BigInt(1), v = BigInt(0);
            while (b != BigInt(0)) {
                BigInt t = a / b;
                a = a - t * b;
                std::swap(a, b);
                u = u - t * v;
                std::swap(u, v);
            }
            if (a != BigInt(1))
                return -1;
            return (u % p + p) % p;
        }

    
    private:
        static BigInt add_abs(const BigInt &a, const BigInt &b) {
            BigInt res;
            res.digits.resize(std::max(a.digits.size(), b.digits.size()) + 1);
            uint64_t carry = 0;
        
            for (size_t i = 0; i < res.digits.size(); ++i) {
                __uint128_t sum = carry;
                if (i < a.digits.size()) sum += a.digits[i];
                if (i < b.digits.size()) sum += b.digits[i];
                res.digits[i] = sum % BASE;
                carry = sum / BASE;
            }
        
            res.remove_leading_zeros();
            return res;
        }
    
        static BigInt sub_abs(const BigInt &a, const BigInt &b) {
            BigInt res = a;
            uint64_t borrow = 0;
        
            for (size_t i = 0; i < res.digits.size(); ++i) {
                int64_t diff = res.digits[i] - borrow - (i < b.digits.size() ? b.digits[i] : 0);
                if (diff < 0) {
                    diff += BASE;
                    borrow = 1;
                } else {
                    borrow = 0;
                }
                res.digits[i] = diff;
            }
        
            res.remove_leading_zeros();
            return res;
        }

        static BigInt naive_mult(const BigInt &a, const BigInt &b) {
            BigInt res;
            res.digits.resize(a.digits.size() + b.digits.size());
            
            for (size_t i = 0; i < a.digits.size(); ++i) {
                uint64_t carry = 0;
                for (size_t j = 0; j < b.digits.size() || carry; ++j) {
                    __uint128_t cur = res.digits[i + j] + 
                                     (__uint128_t)a.digits[i] * (j < b.digits.size() ? b.digits[j] : 0) + carry;
                    res.digits[i + j] = cur % BASE;
                    carry = cur / BASE;
                }
            }
            res.remove_leading_zeros();
            return res;
        }
    
        static BigInt karatsuba_mult(const BigInt &a, const BigInt &b) {
            if (a.digits.size() < 32 || b.digits.size() < 32)
                return naive_mult(a, b);
    
            size_t m = a.digits.size() / 2;
            BigInt a1(std::vector<uint64_t>(a.digits.begin() + m, a.digits.end()));
            BigInt a0(std::vector<uint64_t>(a.digits.begin(), a.digits.begin() + m));
            BigInt b1(std::vector<uint64_t>(b.digits.begin() + m, b.digits.end()));
            BigInt b0(std::vector<uint64_t>(b.digits.begin(), b.digits.begin() + m));
    
            BigInt z2 = karatsuba_mult(a1, b1);
            BigInt z0 = karatsuba_mult(a0, b0);
            BigInt z1 = karatsuba_mult(add_abs(a1, a0), add_abs(b1, b0)) - z2 - z0;
    
            return (z2.shift(2 * m) + z1.shift(m) + z0);
        }
    
        void remove_leading_zeros() {
            while (digits.size() > 1 && digits.back() == 0)
                digits.pop_back();
            if (digits.size() == 1 && digits[0] == 0) is_negative = false;
        }
        
        BigInt shift(size_t m) const {
            BigInt res = *this;
            res.digits.insert(res.digits.begin(), m, 0);
            return res;
        }
    
        BigInt negate() const {
            BigInt res = *this;
            res.is_negative = !is_negative;
            return res;
        }

        static BigInt mul_abs(const BigInt &a, uint64_t b) {
            BigInt res;
            res.digits.resize(a.digits.size() + 1);
        
            __uint128_t carry = 0;
            for (size_t i = 0; i < a.digits.size(); ++i) {
                __uint128_t prod = (__uint128_t)a.digits[i] * b + carry;
                res.digits[i] = prod % BASE;
                carry = prod / BASE;
            }
            res.digits[a.digits.size()] = carry;
            res.remove_leading_zeros();
            return res;
        }
};
    
#endif