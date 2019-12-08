#include <cmath>
#include <climits>
#include "BigInt.h"

BigInt::BigInt() {
    element.resize(50);
    digits = 1;
    status = Status::Plus;
}

BigInt::BigInt(long long n) {
    if(n < 0) {
        status = Status::Minus;
        n *= -1;
    }
    else status = Status::Plus;

    long long l = n;
    l /= 10;
    for(digits = 1; l > 0; digits++) l /= 10;
    element.resize(digits + 50);
    for(int i = 0; i < digits; i++) {
        element[i] = n % 10;
        n /= 10;
    }
}

// 算術演算

BigInt BigInt::operator + () const {
    return *this;
}

BigInt BigInt::operator - () const {
    return *this * (-1);
}

BigInt& BigInt::operator ++ () {
    *this = *this + 1;
    return *this;
}

BigInt BigInt::operator ++ (int n) {
    *this = *this + 1;
    return *this;
}

BigInt& BigInt::operator -- () {
    *this = *this - 1;
    return *this;
}

BigInt BigInt::operator -- (int n) {
    *this = *this - 1;
    return *this;
}

template<typename T>
BigInt operator + (BigInt n, T k) {
    if(n < k) return k + n;

    for(size_t i = 0; k != 0; i++) {
        n.element[i] += k % 10;
        k /= 10;
    }
    for(size_t i = 0; i < n.digits; i++) {
        if(i + 1 == n.element.size())
            n.element.resize(n.element.size() * 2);

        if(i == n.digits - 1 && n.element[i] > 9)
            n.digits++;

        n.element[i + 1] += n.element[i] / 10;
        n.element[i] %= 10;
    }
    return n;
}

template BigInt operator + <char>(BigInt, char);
template BigInt operator + <unsigned char>(BigInt, unsigned char);
template BigInt operator + <short>(BigInt, short);
template BigInt operator + <unsigned short>(BigInt, unsigned short);
template BigInt operator + <int>(BigInt, int);
template BigInt operator + <unsigned int>(BigInt, unsigned int);
template BigInt operator + <long>(BigInt, long);
template BigInt operator + <unsigned long>(BigInt, unsigned long);
template BigInt operator + <long long>(BigInt, long long);

template<typename T>
BigInt operator + (T n, BigInt k) {
    if(n < k) return k + n;

    for(size_t i = 0; n != 0; i++) {
        if(i + 1 == k.element.size())
            k.element.resize(k.element.size() * 2);
        k.element[i] += n % 10;
        n /= 10;
    }
    for(size_t i = 0; i < k.digits; i++) {
        if(i + 1 == k.element.size())
            k.element.resize(k.element.size() * 2);

        if(i == k.digits - 1 && k.element[i] > 9)
            k.digits++;

        k.element[i + 1] += k.element[i] / 10;
        k.element[i] %= 10;
    }
    return k;
}

template BigInt operator + <char>(char, BigInt);
template BigInt operator + <unsigned char>(unsigned char, BigInt);
template BigInt operator + <short>(short, BigInt);
template BigInt operator + <unsigned short>(unsigned short, BigInt);
template BigInt operator + <int>(int, BigInt);
template BigInt operator + <unsigned int>(unsigned int, BigInt);
template BigInt operator + <long>(long, BigInt);
template BigInt operator + <unsigned long>(unsigned long, BigInt);
template BigInt operator + <long long>(long long, BigInt);
template BigInt operator + <unsigned long long>(unsigned long long, BigInt);

BigInt operator + (BigInt n, BigInt k) {
    size_t i;
    if(n < k) return k + n;
    if(k.status == Status::Minus) {
        k.status = Status::Plus;
        return n - k;
    }

    for(i = 0; i < k.digits; i++) {
        n.element[i] += k.element[i];
    }
    for(i = 0; i < n.digits; i++) {
        if(i + 1 == n.element.size())
            n.element.resize(n.element.size() * 2);

        if(i == n.digits - 1 && n.element[i] > 9)
            n.digits++;

        n.element[i + 1] += n.element[i] / 10;
        n.element[i] %= 10;
    }
    return n;
}

template<typename T>
BigInt operator - (BigInt n, T k) {
    if(n < k) return -(k - n);
    if(k < 0) return n + (-k);

    for(size_t i = 0; k != 0; i++) {
        if(n.element[i] < (k % 10)) {
            n.element[i] += 10;
            for(size_t j = 1; ; j++) {
                if(n.element[i + j] != 0) {
                    n.element[i + j]--;
                    break;
                }
                else n.element[i + j] = 9;
            }
        }
        n.element[i] -= k % 10;
        k /= 10;
    }
    
    for(size_t i = n.digits - 1; i > 0; i--) {
        if(n.element[i] == 0) {
            n.digits--;
        }
        else break;
    }
    return n;
}

template BigInt operator - <char>(BigInt, char);
template BigInt operator - <unsigned char>(BigInt, unsigned char);
template BigInt operator - <short>(BigInt, short);
template BigInt operator - <unsigned short>(BigInt, unsigned short);
template BigInt operator - <int>(BigInt, int);
template BigInt operator - <unsigned int>(BigInt, unsigned int);
template BigInt operator - <long>(BigInt, long);
template BigInt operator - <unsigned long>(BigInt, unsigned long);
template BigInt operator - <long long>(BigInt, long long);

template<typename T>
BigInt operator - (T n, BigInt k) {
    return BigInt(n) - k;
}

template BigInt operator - <char>(char, BigInt);
template BigInt operator - <unsigned char>(unsigned char, BigInt);
template BigInt operator - <short>(short, BigInt);
template BigInt operator - <unsigned short>(unsigned short, BigInt);
template BigInt operator - <int>(int, BigInt);
template BigInt operator - <unsigned int>(unsigned int, BigInt);
template BigInt operator - <long>(long, BigInt);
template BigInt operator - <unsigned long>(unsigned long, BigInt);
template BigInt operator - <long long>(long long, BigInt);
template BigInt operator - <unsigned long long>(unsigned long long, BigInt);

BigInt operator - (BigInt n, BigInt k) {
    int i, j, l, tmp, val;
    if(n < k) return -(k - n);
    if(k.status == Status::Minus) {
        k.status = Status::Plus;
        return n + k;
    }
    for(i = 0; i < k.digits; i++) {
        if(n.element[i] < k.element[i]) {
            n.element[i] += 10;
            for(j = 1; ; j++) {
                if(n.element[i + j] != 0) {
                    n.element[i + j]--;
                    break;
                }
                else n.element[i + j] = 9;
            }
        }
        n.element[i] -= k.element[i];
    }
    
    for(i = n.digits - 1; i > 0; i--) {
        if(n.element[i] == 0) {
            n.digits--;
        }
        else break;
    }
    return n;
}

BigInt karatsuba(BigInt n, BigInt k) {
    BigInt n1 = 0, n2 = 0, k1 = 0, k2 = 0;
    BigInt B = 1;
    size_t d = size_t(n.digits / 2);
    for(size_t i = 0; i < d; ++i) {
        n1 += B * n.element[i];
        n2 += B * n.element[i + d];
        if(i < k.digits) k1 += B * k.element[i];
        if(i + d < k.digits) k2 += B * k.element[i + d];
        B *= 10;
    }
    if(2 * d != n.digits) {
        n2 += B * n.element[2 * d];
        if(2 * d < k.digits) k2 += B * k.element[2 * d];
    }
    return n1 * k1 + (n1 * k2 + n2 * k1) * B + n2 * k2 * B * B;
}

template<typename T>
BigInt operator * (BigInt n, T k) {
    return n * BigInt(k);
}

template BigInt operator * <char>(BigInt, char);
template BigInt operator * <unsigned char>(BigInt, unsigned char);
template BigInt operator * <short>(BigInt, short);
template BigInt operator * <unsigned short>(BigInt, unsigned short);
template BigInt operator * <int>(BigInt, int);
template BigInt operator * <unsigned int>(BigInt, unsigned int);
template BigInt operator * <long>(BigInt, long);
template BigInt operator * <unsigned long>(BigInt, unsigned long);
template BigInt operator * <long long>(BigInt, long long);

template<typename T>
BigInt operator * (T n, BigInt k) {
    return BigInt(n) * k;
}

template BigInt operator * <char>(char, BigInt);
template BigInt operator * <unsigned char>(unsigned char, BigInt);
template BigInt operator * <short>(short, BigInt);
template BigInt operator * <unsigned short>(unsigned short, BigInt);
template BigInt operator * <int>(int, BigInt);
template BigInt operator * <unsigned int>(unsigned int, BigInt);
template BigInt operator * <long>(long, BigInt);
template BigInt operator * <unsigned long>(unsigned long, BigInt);
template BigInt operator * <long long>(long long, BigInt);
template BigInt operator * <unsigned long long>(unsigned long long, BigInt);

BigInt operator * (BigInt n, BigInt k) {
    if(n < k) return k * n;
    if(n == 0 || k == 0) return 0;
    // BigInt absn = n.abs(), absk = k.abs();
    // if(absn < int(sqrt(LLONG_MAX)) && absk < int(sqrt(LLONG_MAX))) return BigInt((long long)n * (long long)k);
    // if(absn < int(sqrt(INT_MAX)) && absk < int(sqrt(INT_MAX))) return BigInt(int(n) * int(k));
    size_t i, j;
    BigInt tmp = k;

    if(n.status == k.status)
        k.status = Status::Plus;
    else k.status = Status::Minus;

    if(k.digits == k.element.size()) k.element.resize(k.element.size() * 2);

    if(n.digits > 2 && k.digits > 2) {
        while(k.element[0] == 0) {
            k.element.erase(k.element.begin());
            n.element.insert(n.element.begin(), 0);
            k.digits--;
            n.digits++;
        }
        if(k.digits > 2) return karatsuba(n, k);
        else return n * k;
    }

    for(i = 0; i < k.digits; i++)
        k.element[i] = 0;
    
    for(i = 0; i < tmp.digits; i++) {
        for(j = 0; j < n.digits; j++) {
            if(i + j == k.element.size())
                k.element.resize(k.element.size() * 2);
            
            if(i + j + 1 > k.digits)
                k.digits = i + j + 1;

            k.element[i + j] += tmp.element[i] * n.element[j];
        }
    }

    for(i = 0; i < k.digits; i++) {
        if(i + 1 == k.element.size())
            k.element.resize(k.element.size() * 2);
        
        if(i == k.digits - 1 && k.element[i] > 9)
            k.digits++;
        
        k.element[i + 1] += k.element[i] / 10;
        k.element[i] %= 10;
    }
    return k;
}

template<typename T>
BigInt operator / (BigInt n, T k) {
    return n / BigInt(k);
}

template BigInt operator / <char>(BigInt, char);
template BigInt operator / <unsigned char>(BigInt, unsigned char);
template BigInt operator / <short>(BigInt, short);
template BigInt operator / <unsigned short>(BigInt, unsigned short);
template BigInt operator / <int>(BigInt, int);
template BigInt operator / <unsigned int>(BigInt, unsigned int);
template BigInt operator / <long>(BigInt, long);
template BigInt operator / <unsigned long>(BigInt, unsigned long);
template BigInt operator / <long long>(BigInt, long long);

template<typename T>
BigInt operator / (T n, BigInt k) {
    return BigInt(n) / k;
}

template BigInt operator / <char>(char, BigInt);
template BigInt operator / <unsigned char>(unsigned char, BigInt);
template BigInt operator / <short>(short, BigInt);
template BigInt operator / <unsigned short>(unsigned short, BigInt);
template BigInt operator / <int>(int, BigInt);
template BigInt operator / <unsigned int>(unsigned int, BigInt);
template BigInt operator / <long>(long, BigInt);
template BigInt operator / <unsigned long>(unsigned long, BigInt);
template BigInt operator / <long long>(long long, BigInt);
template BigInt operator / <unsigned long long>(unsigned long long, BigInt);

BigInt operator / (BigInt n, BigInt k) {
    if(k == 0) throw "Divide by Zero";
    BigInt c = 1, q = 0;
    BigInt absn = n.abs(), absk = k.abs();
    BigInt save_k = absk;
    while(absn >= absk) {
        if(absn >= 2 * absk) {
            absk *= 2;
            c *= 2;
        }
        else {
            q += c;
            c = 1;
            absn -= absk;
            absk = save_k;
        }
    }
    return q * int(n.status) * int(k.status);
}

template<typename T>
BigInt operator % (BigInt n, T k) {
    return n % BigInt(k);
}

template BigInt operator % <char>(BigInt, char);
template BigInt operator % <unsigned char>(BigInt, unsigned char);
template BigInt operator % <short>(BigInt, short);
template BigInt operator % <unsigned short>(BigInt, unsigned short);
template BigInt operator % <int>(BigInt, int);
template BigInt operator % <unsigned int>(BigInt, unsigned int);
template BigInt operator % <long>(BigInt, long);
template BigInt operator % <unsigned long>(BigInt, unsigned long);
template BigInt operator % <long long>(BigInt, long long);

template<typename T>
BigInt operator % (T n, BigInt k) {
    return BigInt(n) % k;
}

template BigInt operator % <char>(char, BigInt);
template BigInt operator % <unsigned char>(unsigned char, BigInt);
template BigInt operator % <short>(short, BigInt);
template BigInt operator % <unsigned short>(unsigned short, BigInt);
template BigInt operator % <int>(int, BigInt);
template BigInt operator % <unsigned int>(unsigned int, BigInt);
template BigInt operator % <long>(long, BigInt);
template BigInt operator % <unsigned long>(unsigned long, BigInt);
template BigInt operator % <long long>(long long, BigInt);
template BigInt operator % <unsigned long long>(unsigned long long, BigInt);

BigInt operator % (BigInt n, BigInt k) {
    return n - k * (n / k);
}

BigInt& BigInt::operator ~ () {
    std::vector<size_t> b = this->to_binary();
    for(int i = 0; i < b.size(); i++) {
        b[i] ^= 1;
    }
    *this = binary_to_i(b);
    return *this;
}

BigInt operator & (BigInt n, BigInt k) {
    if(n < k) return k & n;
    std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
    int i;
    for(i = 0; i < b_k.size(); i++) {
        b_k[i] &= b_n[i];
    }
    return binary_to_i(b_k);
}

BigInt operator | (BigInt n, BigInt k) {
    if(n < k) return k | n;
    std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
    int i;
    for(i = 0; i < b_k.size(); i++) {
        b_n[i] |= b_k[i];
    }
    return binary_to_i(b_n);
}

BigInt operator ^ (BigInt n, BigInt k) {
    if(n < k) return k ^ n;
    std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
    int i;
    for(i = 0; i < b_k.size(); i++) {
        b_n[i] ^= b_k[i];
    }
    return binary_to_i(b_n);
}

template<typename T>
BigInt BigInt::operator << (T n) {
    std::vector<size_t> b = this->to_binary();
    for(size_t i = 0; i < n; ++i) {
        b.insert(b.begin(), 0);
    }
    return binary_to_i(b);
}

template BigInt BigInt::operator << <int>(int);

template<typename T>
BigInt BigInt::operator >> (T n) {
    std::vector<size_t> b = this->to_binary();
    b.erase(b.begin(), b.begin() + n);
    return binary_to_i(b);
}

template BigInt BigInt::operator >> <int>(int);

// 代入演算

BigInt& BigInt::operator = (long long n) {
    this->digits = sizeof(n) / sizeof(long long);
    for(int i = 0; n > 0; i++) {
        if(i == this->element.size())
            this->element.resize(this->element.size() * 2);
        
        this->element[i] = n % 10;
        n /= 10;
    }
    return *this;
}

BigInt& BigInt::operator += (const BigInt n) {
    *this = *this + n;
    return *this;
}

BigInt& BigInt::operator -= (const BigInt n) {
    *this = *this - n;
    return *this;
}

BigInt& BigInt::operator *= (const BigInt n) {
    *this = *this * n;
    return *this;
}

BigInt& BigInt::operator /= (const BigInt n) {
    *this = *this / n;
    return *this;
}

BigInt& BigInt::operator %= (const BigInt n) {
    *this = *this % n;
    return *this;
}

BigInt& BigInt::operator &= (const BigInt n) {
    *this = *this & n;
    return *this;
}

BigInt& BigInt::operator |= (const BigInt n) {
    *this = *this | n;
    return *this;
}

BigInt& BigInt::operator ^= (const BigInt n) {
    *this = *this ^ n;
    return *this;
}

BigInt& BigInt::operator <<= (const size_t n) {
    *this = *this << n;
    return *this;
}

BigInt& BigInt::operator >>= (const size_t n) {
    *this = *this >> n;
    return *this;
}

// 比較演算

template<typename T>
bool operator > (BigInt n, T k) {
    return n > BigInt(k);
}

template bool operator > <char>(BigInt, char);
template bool operator > <unsigned char>(BigInt, unsigned char);
template bool operator > <short>(BigInt, short);
template bool operator > <unsigned short>(BigInt, unsigned short);
template bool operator > <int>(BigInt, int);
template bool operator > <unsigned int>(BigInt, unsigned int);
template bool operator > <long>(BigInt, long);
template bool operator > <unsigned long>(BigInt, unsigned long);
template bool operator > <long long>(BigInt, long long);

template<typename T>
bool operator > (T n, BigInt k) {
    return BigInt(n) > k;
}

template bool operator > <char>(char, BigInt);
template bool operator > <unsigned char>(unsigned char, BigInt);
template bool operator > <short>(short, BigInt);
template bool operator > <unsigned short>(unsigned short, BigInt);
template bool operator > <int>(int, BigInt);
template bool operator > <unsigned int>(unsigned int, BigInt);
template bool operator > <long>(long, BigInt);
template bool operator > <unsigned long>(unsigned long, BigInt);
template bool operator > <long long>(long long, BigInt);
template bool operator > <unsigned long long>(unsigned long long, BigInt);

bool operator > (BigInt n, BigInt k) {
    if(n.status == Status::Plus && k.status == Status::Minus) return true;
    if(n.status == Status::Minus && k.status == Status::Plus) return false;

    bool flag = true;
    if(n.status == Status::Minus) flag = false;

    if(n.digits > k.digits) return flag;
    if(n.digits < k.digits) return !flag;

    int i;
    for(i = n.digits - 1; i >= 0; i--) {
        if(n.element[i] > k.element[i]) return flag;
        if(n.element[i] < k.element[i]) return !flag;
    }
    return false;
}

template<typename T>
bool operator >= (BigInt n, T k) {
    return n >= BigInt(k);
}

template bool operator >= <char>(BigInt, char);
template bool operator >= <unsigned char>(BigInt, unsigned char);
template bool operator >= <short>(BigInt, short);
template bool operator >= <unsigned short>(BigInt, unsigned short);
template bool operator >= <int>(BigInt, int);
template bool operator >= <unsigned int>(BigInt, unsigned int);
template bool operator >= <long>(BigInt, long);
template bool operator >= <unsigned long>(BigInt, unsigned long);
template bool operator >= <long long>(BigInt, long long);

template<typename T>
bool operator >= (T n, BigInt k) {
    return BigInt(n) >= k;
}

template bool operator >= <char>(char, BigInt);
template bool operator >= <unsigned char>(unsigned char, BigInt);
template bool operator >= <short>(short, BigInt);
template bool operator >= <unsigned short>(unsigned short, BigInt);
template bool operator >= <int>(int, BigInt);
template bool operator >= <unsigned int>(unsigned int, BigInt);
template bool operator >= <long>(long, BigInt);
template bool operator >= <unsigned long>(unsigned long, BigInt);
template bool operator >= <long long>(long long, BigInt);
template bool operator >= <unsigned long long>(unsigned long long, BigInt);

bool operator >= (BigInt n, BigInt k) {
    if(n.status == Status::Plus && k.status == Status::Minus) return true;
    if(n.status == Status::Minus && k.status == Status::Plus) return false;

    bool flag = true;
    if(n.status == Status::Minus) flag = false;

    if(n.digits > k.digits) return flag;
    if(n.digits < k.digits) return !flag;

    int i;
    for(i = n.digits - 1; i >= 0; i--) {
        if(n.element[i] > k.element[i]) return flag;
        if(n.element[i] < k.element[i]) return !flag;
    }
    return true;
}

template<typename T>
bool operator < (BigInt n, T k) {
    return n < BigInt(k);
}

template bool operator < <char>(BigInt, char);
template bool operator < <unsigned char>(BigInt, unsigned char);
template bool operator < <short>(BigInt, short);
template bool operator < <unsigned short>(BigInt, unsigned short);
template bool operator < <int>(BigInt, int);
template bool operator < <unsigned int>(BigInt, unsigned int);
template bool operator < <long>(BigInt, long);
template bool operator < <unsigned long>(BigInt, unsigned long);
template bool operator < <long long>(BigInt, long long);

template<typename T>
bool operator < (T n, BigInt k) {
    return BigInt(n) < k;
}

template bool operator < <char>(char, BigInt);
template bool operator < <unsigned char>(unsigned char, BigInt);
template bool operator < <short>(short, BigInt);
template bool operator < <unsigned short>(unsigned short, BigInt);
template bool operator < <int>(int, BigInt);
template bool operator < <unsigned int>(unsigned int, BigInt);
template bool operator < <long>(long, BigInt);
template bool operator < <unsigned long>(unsigned long, BigInt);
template bool operator < <long long>(long long, BigInt);
template bool operator < <unsigned long long>(unsigned long long, BigInt);

bool operator < (BigInt n, BigInt k) {
    if(n.status == Status::Plus && k.status == Status::Minus) return false;
    if(n.status == Status::Minus && k.status == Status::Plus) return true;

    bool flag = true;
    if(n.status == Status::Minus) flag = false;

    if(n.digits < k.digits) return flag;
    if(n.digits > k.digits) return !flag;

    int i;
    for(i = n.digits - 1; i >= 0; i--) {
        if(n.element[i] < k.element[i]) return flag;
        if(n.element[i] > k.element[i]) return !flag;
    }
    return false;
}

template<typename T>
bool operator <= (BigInt n, T k) {
    return n <= BigInt(k);
}

template bool operator <= <char>(BigInt, char);
template bool operator <= <unsigned char>(BigInt, unsigned char);
template bool operator <= <short>(BigInt, short);
template bool operator <= <unsigned short>(BigInt, unsigned short);
template bool operator <= <int>(BigInt, int);
template bool operator <= <unsigned int>(BigInt, unsigned int);
template bool operator <= <long>(BigInt, long);
template bool operator <= <unsigned long>(BigInt, unsigned long);
template bool operator <= <long long>(BigInt, long long);

template<typename T>
bool operator <= (T n, BigInt k) {
    return BigInt(n) <= k;
}

template bool operator <= <char>(char, BigInt);
template bool operator <= <unsigned char>(unsigned char, BigInt);
template bool operator <= <short>(short, BigInt);
template bool operator <= <unsigned short>(unsigned short, BigInt);
template bool operator <= <int>(int, BigInt);
template bool operator <= <unsigned int>(unsigned int, BigInt);
template bool operator <= <long>(long, BigInt);
template bool operator <= <unsigned long>(unsigned long, BigInt);
template bool operator <= <long long>(long long, BigInt);
template bool operator <= <unsigned long long>(unsigned long long, BigInt);

bool operator <= (BigInt n, BigInt k) {
    if(n.status == Status::Plus && k.status == Status::Minus) return false;
    if(n.status == Status::Minus && k.status == Status::Plus) return true;

    bool flag = true;
    if(n.status == Status::Minus) flag = false;

    if(n.digits < k.digits) return flag;
    if(n.digits > k.digits) return !flag;

    int i;
    for(i = n.digits - 1; i >= 0; i--) {
        if(n.element[i] < k.element[i]) return flag;
        if(n.element[i] > k.element[i]) return !flag;
    }
    return true;
}

template<typename T>
bool operator == (BigInt n, T k) {
    return n == BigInt(k);
}

template bool operator == <char>(BigInt, char);
template bool operator == <unsigned char>(BigInt, unsigned char);
template bool operator == <short>(BigInt, short);
template bool operator == <unsigned short>(BigInt, unsigned short);
template bool operator == <int>(BigInt, int);
template bool operator == <unsigned int>(BigInt, unsigned int);
template bool operator == <long>(BigInt, long);
template bool operator == <unsigned long>(BigInt, unsigned long);
template bool operator == <long long>(BigInt, long long);

template<typename T>
bool operator == (T n, BigInt k) {
    return BigInt(n) == k;
}

template bool operator == <char>(char, BigInt);
template bool operator == <unsigned char>(unsigned char, BigInt);
template bool operator == <short>(short, BigInt);
template bool operator == <unsigned short>(unsigned short, BigInt);
template bool operator == <int>(int, BigInt);
template bool operator == <unsigned int>(unsigned int, BigInt);
template bool operator == <long>(long, BigInt);
template bool operator == <unsigned long>(unsigned long, BigInt);
template bool operator == <long long>(long long, BigInt);
template bool operator == <unsigned long long>(unsigned long long, BigInt);

bool operator == (BigInt n, BigInt k) {
    if(n.status != k.status) return false;
    if(n.digits != k.digits) return false;
    for(int i = 0; i < n.digits; i++) {
        if(n.element[i] != k.element[i]) return false;
    }
    return true;
}

template<typename T>
bool operator != (BigInt n, T k) {
    return n != BigInt(k);
}

template bool operator != <char>(BigInt, char);
template bool operator != <unsigned char>(BigInt, unsigned char);
template bool operator != <short>(BigInt, short);
template bool operator != <unsigned short>(BigInt, unsigned short);
template bool operator != <int>(BigInt, int);
template bool operator != <unsigned int>(BigInt, unsigned int);
template bool operator != <long>(BigInt, long);
template bool operator != <unsigned long>(BigInt, unsigned long);
template bool operator != <long long>(BigInt, long long);

template<typename T>
bool operator != (T n, BigInt k) {
    return BigInt(n) != k;
}

template bool operator != <char>(char, BigInt);
template bool operator != <unsigned char>(unsigned char, BigInt);
template bool operator != <short>(short, BigInt);
template bool operator != <unsigned short>(unsigned short, BigInt);
template bool operator != <int>(int, BigInt);
template bool operator != <unsigned int>(unsigned int, BigInt);
template bool operator != <long>(long, BigInt);
template bool operator != <unsigned long>(unsigned long, BigInt);
template bool operator != <long long>(long long, BigInt);
template bool operator != <unsigned long long>(unsigned long long, BigInt);

bool operator != (BigInt n, BigInt k) {
    if(n.status != k.status) return true;
    if(n.digits != k.digits) return true;
    for(int i = 0; i < n.digits; i++) {
        if(n.element[i] != k.element[i]) return true;
    }
    return false;
}

// インデックス

BigInt BigInt::operator [] (const BigInt n) {
    auto itr = this->element.begin();
    for(BigInt i = 0; i < n; i++) {
        itr++;
    }
    return *itr;
}

// 入力・出力

std::istream& operator >> (std::istream& in, BigInt& n) {
    long long k;
    in >> k;
    if(k < 0) {
        n.status = Status::Minus;
        k *= -1;
    }
    for(int i = 0; k > 0; i++) {
        if(i == n.element.size())
            n.element.resize(n.element.size() * 2);
        
        n.element[i] = k % 10;
        k /= 10;
        n.digits = i + 1;
    }
    return in;
}

std::ostream& operator << (std::ostream& out, BigInt n) {
    if(n.status == Status::Minus) out << '-';
    for(int i = n.digits - 1; i >= 0; i--) {
        out << n.element[i];
    }
    return out;
}

// キャスト

BigInt::operator char() const {
    if(*this > CHAR_MAX || *this < CHAR_MIN) return 0;
    char c = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        c += a * this->element[i];
        a *= 10;
    }
    return c * char(this->status);
}

BigInt::operator unsigned char() const {
    if(*this < 0 || *this > UCHAR_MAX) return 0;
    unsigned char uc = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        uc += a * this->element[i];
        a *= 10;
    }
    return uc;
}

BigInt::operator short() const {
    if(*this > SHRT_MAX || *this < SHRT_MIN) return 0;
    short s = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        s += a * this->element[i];
        a *= 10;
    }
    return s * short(this->status);
}

BigInt::operator unsigned short() const {
    if(*this < 0 || *this > USHRT_MAX) return 0;
    unsigned short us = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        us += a * this->element[i];
        a *= 10;
    }
    return us;
}

BigInt::operator int() const {
    if(*this > INT_MAX || *this < INT_MIN) return 0;
    int i = 0, a = 1;
    for(size_t j = 0; j < this->digits; ++j) {
        i += a * this->element[j];
        a *= 10;
    }
    return i * int(this->status);
}

BigInt::operator unsigned int() const {
    if(*this < 0 || *this > UINT_MAX) return 0;
    unsigned int ui = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        ui += a * this->element[i];
        a *= 10;
    }
    return ui;
}

BigInt::operator long() const {
    if(*this > LONG_MAX || *this < LONG_MIN) return 0;
    long l = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        l += a * this->element[i];
        a *= 10;
    }
    return l * int(this->status);
}

BigInt::operator unsigned long() const {
    if(*this < 0 || *this > ULONG_MAX) return 0;
    unsigned long l = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        l += a * this->element[i];
        a *= 10;
    }
    return l;
}

BigInt::operator long long() const {
    if(*this > LLONG_MAX || *this < LLONG_MIN) return 0;
    long ll = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        ll += a * this->element[i];
        a *= 10;
    }
    return ll * int(this->status);
}

BigInt::operator unsigned long long() const {
    if(*this < 0 || *this > ULLONG_MAX) return 0;
    unsigned long ll = 0, a = 1;
    for(size_t i = 0; i < this->digits; ++i) {
        ll += a * this->element[i];
        a *= 10;
    }
    return ll;
}

// BigInt::operator float() const {
//     if(*this > LC_MAX || *this < LC_MIN) return 0.0;
//     float f = 0.0;
//     int a = 1;
//     for(size_t i = 0; i < this->digits; ++i) {
//         f += a * this->element[i];
//         a *= 10;
//     }
//     return f;
// }

// BigInt::operator double() const {
//     if(*this > TMP_MAX || *this < TMP_MAX) return 0.0;
//     double d = 0.0;
//     int a = 1;
//     for(size_t i = 0; i < this->digits; ++i) {
//         d += a * this->element[i];
//         a *= 10;
//     }
//     return d;
// }

// その他

size_t BigInt::size() {
    return digits;
}

std::vector<size_t> BigInt::to_binary() {
    std::vector<size_t> binary;
    BigInt n = *this;
    while(n > 0) {
        if(n % 2 == 0) binary.push_back(0);
        else binary.push_back(1);
        n /= 2;
    }
    return binary;
}

BigInt BigInt::abs() {
    return *this * int(this->status);
}

BigInt binary_to_i(std::vector<size_t> b) {
    int i, j, val;
    BigInt n;
    for(i = 0; i < b.size(); i++) {
        val = b[i];
        for(j = 0; j < i; j++) {
            val *= 2;
        }
        n += val;
    }
    return n;
}