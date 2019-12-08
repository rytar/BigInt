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
    return *this + 1;
}

BigInt& BigInt::operator -- () {
    *this = *this - 1;
    return *this;
}

BigInt BigInt::operator -- (int n) {
    return *this - 1;
}

BigInt operator + (BigInt n, BigInt k) {
    int i;
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

BigInt operator * (BigInt n, BigInt k) {
    int i, j;
    BigInt tmp = k;
    if(n == 0 || k == 0) return 0;

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

BigInt BigInt::operator << (size_t n) {
    std::vector<size_t> b = this->to_binary();
    for(size_t i = 0; i < n; ++i) {
        b.insert(b.begin(), 0);
    }
    return binary_to_i(b);
}

BigInt BigInt::operator >> (size_t n) {
    std::vector<size_t> b = this->to_binary();
    b.erase(b.begin(), b.begin() + n);
    return binary_to_i(b);
}

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


bool operator == (BigInt n, BigInt k) {
    if(n.status != k.status) return false;
    if(n.digits != k.digits) return false;
    for(int i = 0; i < n.digits; i++) {
        if(n.element[i] != k.element[i]) return false;
    }
    return true;
}

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

// BigInt::operator char() const {
//     char c = 0;
//     size_t a = 1;
//     for(size_t i = 0; i < this->digits; i++) {
//         c += a * this->element[i];
//         a *= 10;
//     }
//     return c;
// }

// BigInt::operator unsigned char() const {
//     unsigned char uc = 0;
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




