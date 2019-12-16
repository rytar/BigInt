#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>

namespace Rytar {

    enum class Status {
        Plus = 1, Minus = -1
    };

    class BigInt {
    private:
        std::vector<size_t> element;
        size_t digits;
        Status status;

    public:
        BigInt() {
            element.resize(50);
            digits = 1;
            status = Status::Plus;
        }

        template<typename T>
        BigInt(T n) {
            if(n < 0) {
                status = Status::Minus;
                n *= -1;
            }
            else status = Status::Plus;

            digits = 0;
            while(n > 0) {
                element.push_back(n % 10);
                n /= 10;
                digits++;
            }
            element.resize(digits + 50);
        }

        BigInt(const char* c) {
            auto s = std::string(c);
            digits = s.size();
            element.resize(digits + 50);
            status = Status::Plus;
            for(size_t i = 0; i < digits + 0; ++i) {
                if(s[digits - i - 1] == '-') status = Status::Minus;
                else element[i] = s[digits - i - 1] - '0';
            }
        }
        
        BigInt(const std::string s) {
            digits = s.size();
            element.resize(digits + 50);
            status = Status::Plus;
            for(size_t i = 0; i < digits + 0; ++i) {
                if(s[digits - i - 1] == '-') status = Status::Minus;
                else element[i] = s[digits - i - 1] - '0';
            }
        }

        // 算術演算
        BigInt operator + () const {
            return *this;
        }

        BigInt operator - () const {
            return *this * (-1);
        }

        BigInt& operator ++ () {
            *this = *this + 1;
            return *this;
        }

        BigInt operator ++ (int) {
            *this = *this + 1;
            return *this;
        }

        BigInt& operator -- () {
            *this = *this - 1;
            return *this;
        }

        BigInt operator -- (int) {
            *this = *this - 1;
            return *this;
        }

        template<typename T>
        friend BigInt operator + (BigInt, T);
        template<typename T>
        friend BigInt operator + (T, BigInt);
        friend BigInt operator + (BigInt, BigInt);
        template<typename T>
        friend BigInt operator - (BigInt, T);
        template<typename T>
        friend BigInt operator - (T, BigInt);
        friend BigInt operator - (BigInt, BigInt);
        
        static BigInt karatsuba(BigInt n, BigInt k) {
            BigInt n1 = 0, n2 = 0, k1 = 0, k2 = 0;
            BigInt B = 1;
            size_t d = n.digits / 2;
            for(size_t i = 0; i < d; ++i) {
                n1 += B * n.element[i];
                n2 += B * n.element[i + d];
                if(i < k.digits) k1 += B * k.element[i];
                if(i + d < k.digits) k2 += B * k.element[i + d];
                B *= 10;
            }
            if(2 * d < n.digits) n2 += B * n.element[2 * d];
            if(2 * d < k.digits) k2 += B * k.element[2 * d];
            return n1 * k1 + (n1 * k2 + n2 * k1) * B + n2 * k2 * B * B;
        }
        
        template<typename T>
        friend BigInt operator * (BigInt, T);
        template<typename T>
        friend BigInt operator * (T, BigInt);
        friend BigInt operator * (BigInt, BigInt);
        
        static BigInt karatsuba_for_div(BigInt n, BigInt k) {
            BigInt n1 = 0, n2 = 0;
            BigInt B = 1;
            size_t d = n.digits / 2;
            for(size_t i = 0; i < d; ++i) {
                n1 += B * n.element[i];
                n2 += B * n.element[i + d];
                B *= 10;
            }
            if(2 * d < n.digits) n2 += B * n.element[2 * d];
            BigInt q1 = n1 / k, q2 = n2 / k;
            return q1 + q2 * B + ((n1 - k * q1) + (n2 - k * q2) * B) / k;
        }
        
        template<typename T>
        friend BigInt operator / (BigInt, T);
        template<typename T>
        friend BigInt operator / (T, BigInt);
        friend BigInt operator / (BigInt, BigInt);
        template<typename T>
        friend BigInt operator % (BigInt, T);
        template<typename T>
        friend BigInt operator % (T, BigInt);
        friend BigInt operator % (BigInt, BigInt);


        // 関係演算
        template<typename T>
        friend bool operator > (BigInt, T);
        template<typename T>
        friend bool operator > (T, BigInt);
        friend bool operator > (BigInt, BigInt);
        template<typename T>
        friend bool operator >= (BigInt, T);
        template<typename T>
        friend bool operator >= (T, BigInt);
        friend bool operator >= (BigInt, BigInt);
        template<typename T>
        friend bool operator < (BigInt, T);
        template<typename T>
        friend bool operator < (T, BigInt);
        friend bool operator < (BigInt, BigInt);
        template<typename T>
        friend bool operator <= (BigInt, T);
        template<typename T>
        friend bool operator <= (T, BigInt);
        friend bool operator <= (BigInt, BigInt);
        template<typename T>
        friend bool operator == (BigInt, T);
        template<typename T>
        friend bool operator == (T, BigInt);
        friend bool operator == (BigInt, BigInt);
        template<typename T>
        friend bool operator != (BigInt, T);
        template<typename T>
        friend bool operator != (T, BigInt);
        friend bool operator != (BigInt, BigInt);


        // 論理演算
        friend BigInt operator ~ (BigInt);
        template<typename T>
        friend BigInt operator & (BigInt, T);
        template<typename T>
        friend BigInt operator & (T, BigInt);
        friend BigInt operator & (BigInt, BigInt);
        template<typename T>
        friend BigInt operator | (BigInt, T);
        template<typename T>
        friend BigInt operator | (T, BigInt);
        friend BigInt operator | (BigInt, BigInt);
        template<typename T>
        friend BigInt operator ^ (BigInt, T);
        template<typename T>
        friend BigInt operator ^ (T, BigInt);
        friend BigInt operator ^ (BigInt, BigInt);
        BigInt operator << (size_t);
        BigInt operator >> (size_t);


        // 代入演算
        template<typename T>
        BigInt& operator = (const T n) {
            *this = BigInt(n);
            return *this;
        }

        BigInt& operator = (const BigInt n) {
            if(this->digits < n.digits) this->element.resize(n.digits + 50);
            this->status = n.status;
            this->digits = n.digits;
            for(size_t i = 0; i < n.digits; i++)
                this->element[i] = n.element[i];
            return *this;
        }

        template<typename T>
        BigInt& operator += (const T n) {
            *this = *this + BigInt(n);
            return *this;
        }

        BigInt& operator += (const BigInt n) {
            *this = *this + n;
            return *this;
        }

        template<typename T>
        BigInt& operator -= (const T n) {
            *this = *this - BigInt(n);
            return *this;
        }

        BigInt& operator -= (const BigInt n) {
            *this = *this - n;
            return *this;
        }

        template<typename T>
        BigInt& operator *= (const T n) {
            *this = *this * BigInt(n);
            return *this;
        }

        BigInt& operator *= (const BigInt n) {
            *this = *this * n;
            return *this;
        }

        template<typename T>
        BigInt& operator /= (const T n) {
            *this = *this / BigInt(n);
            return *this;
        }

        BigInt& operator /= (const BigInt n) {
            *this = *this / n;
            return *this;
        }

        template<typename T>
        BigInt& operator %= (const T n) {
            *this = *this % BigInt(n);
            return *this;
        }

        BigInt& operator %= (const BigInt n) {
            *this = *this % n;
            return *this;
        }

        template<typename T>
        BigInt& operator &= (const T n) {
            *this = *this & BigInt(n);
            return *this;
        }

        BigInt& operator &= (const BigInt n) {
            *this = *this & n;
            return *this;
        }
        
        template<typename T>
        BigInt& operator |= (const T n) {
            *this = *this | BigInt(n);
            return *this;
        }

        BigInt& operator |= (const BigInt n) {
            *this = *this | n;
            return *this;
        }

        template<typename T>
        BigInt& operator ^= (const T n) {
            *this = *this ^ BigInt(n);
            return *this;
        }

        BigInt& operator ^= (const BigInt n) {
            *this = *this ^ n;
            return *this;
        }

        BigInt& operator <<= (const size_t n) {
            *this = *this << n;
            return *this;
        }

        BigInt& operator >>= (const size_t n) {
            *this = *this >> n;
            return *this;
        }

        // インデックス
        size_t operator [] (const size_t n) {
            auto itr = this->element.begin();
            itr += n;
            return *itr;
        }

        // 入力・出力
        friend std::istream& operator >> (std::istream&, BigInt&);
        friend std::ostream& operator << (std::ostream&, const BigInt);

        // キャスト
        char to_char() const {
            if(*this > CHAR_MAX || *this < CHAR_MIN) return 0;
            char c = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                c += a * this->element[i];
                a *= 10;
            }
            return c * (char)this->status;
        }

        unsigned char to_unsigned_char() const {
            if(*this < 0 || *this > UCHAR_MAX) return 0;
            unsigned char uc = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                uc += a * this->element[i];
                a *= 10;
            }
            return uc;
        }

        short to_short() const {
            if(*this > SHRT_MAX || *this < SHRT_MIN) return 0;
            short s = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                s += a * this->element[i];
                a *= 10;
            }
            return s * (short)this->status;
        }

        unsigned short to_unsigned_short() const {
            if(*this < 0 || *this > USHRT_MAX) return 0;
            unsigned short us = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                us += a * this->element[i];
                a *= 10;
            }
            return us;
        }

        int to_int() const {
            if(*this > INT_MAX || *this < INT_MIN) return 0;
            int i = 0, a = 1;
            for(size_t j = 0; j < this->digits; ++j) {
                i += a * this->element[j];
                a *= 10;
            }
            return i * (int)this->status;
        }

        unsigned int to_unsigned_int() const {
            if(*this < 0 || *this > UINT_MAX) return 0;
            unsigned int ui = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                ui += a * this->element[i];
                a *= 10;
            }
            return ui;
        }

        long to_long() const {
            if(*this > LONG_MAX || *this < LONG_MIN) return 0;
            long l = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                l += a * this->element[i];
                a *= 10;
            }
            return l * (long)this->status;
        }

        unsigned long to_unsigned_long() const {
            if(*this < 0 || *this > ULONG_MAX) return 0;
            unsigned long ul = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                ul += a * this->element[i];
                a *= 10;
            }
            return ul;
        }

        long long to_long_long() const {
            if(*this > LLONG_MAX || *this < LLONG_MIN) return 0;
            long long ll = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                ll += a * this->element[i];
                a *= 10;
            }
            return ll * (long long)this->status;
        }

        unsigned long long to_unsigned_long_long() const {
            if(*this < 0 || *this > ULLONG_MAX) return 0;
            unsigned long long ull = 0, a = 1;
            for(size_t i = 0; i < this->digits; ++i) {
                ull += a * this->element[i];
                a *= 10;
            }
            return ull;
        }

        // その他
        size_t size() {
            return digits;
        }

        BigInt abs() {
            return *this * int(this->status);
        }

        std::vector<size_t> to_binary() {
            std::vector<size_t> binary;
            BigInt n = *this;
            while(n > 0) {
                binary.push_back((n % 2).to_unsigned_int());
                n /= 2;
            }
            return binary;
        }
    };

    BigInt binary_to_i(std::vector<size_t> b) {
        size_t i, j, val;
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

    // 算術演算
    template<typename T>
    BigInt operator + (BigInt n, T k) {
        return n + BigInt(k);
    }

    template<typename T>
    BigInt operator + (T n, BigInt k) {
        return BigInt(n) + k;
    }

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
        return n - BigInt(k);
    }

    template<typename T>
    BigInt operator - (T n, BigInt k) {
        return BigInt(n) - k;
    }

    BigInt operator - (BigInt n, BigInt k) {
        size_t i, j, l, tmp, val;
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

    template<typename T>
    BigInt operator * (BigInt n, T k) {
        return n * BigInt(k);
    }

    template<typename T>
    BigInt operator * (T n, BigInt k) {
        return BigInt(n) * k;
    }

    BigInt operator * (BigInt n, BigInt k) {
        if(n < k) return k * n;
        if(n == 0 || k == 0) return 0;
        size_t i, j;
        BigInt tmp = k;

        if(n.status == k.status)
            k.status = Status::Plus;
        else k.status = Status::Minus;

        k.element.resize(k.element.size() + 10);

        if(n.digits > 2 && k.digits > 2) {
            while(k.element[0] == 0) {
                k.element.erase(k.element.begin());
                n.element.insert(n.element.begin(), 0);
                k.digits--;
                n.digits++;
            }
            if(k.digits > 2) return BigInt::karatsuba(n, k);
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

    template<typename T>
    BigInt operator / (T n, BigInt k) {
        return BigInt(n) / k;
    }

    BigInt operator / (BigInt n, BigInt k) {
        if(k == 0) throw "Divide by Zero";
        if(k == 1 || k == -1) return n * int(k.status);
        if(n < k) return 0;
        BigInt absn = n.abs(), absk = k.abs(), save_k = k.abs();
        if(absn < CHAR_MAX) return n.to_char() / k.to_char();
        if(absn < SHRT_MAX) return n.to_short() / k.to_short();
        if(absn < INT_MAX) return n.to_int() / k.to_int();
        if(absn < LONG_MAX) return n.to_long() / k.to_long();
        if(absn < LLONG_MAX) return n.to_long_long() / k.to_long_long();
        if(n.digits > 2 * k.digits + 2) return BigInt::karatsuba_for_div(absn, absk) * int(n.status) * int(k.status);
        BigInt count = 0;
        while(absn >= absk) {
            absk += save_k;
            count++;
        }
        return count * int(n.status) * int(k.status);
    }

    template<typename T>
    BigInt operator % (BigInt n, T k) {
        return n % BigInt(k);
    }

    template<typename T>
    BigInt operator % (T n, BigInt k) {
        return BigInt(n) % k;
    }

    BigInt operator % (BigInt n, BigInt k) {
        return n - k * (n / k);
    }

    // 関係演算
    template<typename T>
    bool operator > (BigInt n, T k) {
        return n > BigInt(k);
    }

    template<typename T>
    bool operator > (T n, BigInt k) {
        return BigInt(n) > k;
    }

    bool operator > (BigInt n, BigInt k) {
        if(n.status == Status::Plus && k.status == Status::Minus) return true;
        if(n.status == Status::Minus && k.status == Status::Plus) return false;

        bool flag = true;
        if(n.status == Status::Minus) flag = false;

        if(n.digits > k.digits) return flag;
        if(n.digits < k.digits) return !flag;

        size_t i;
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

    template<typename T>
    bool operator >= (T n, BigInt k) {
        return BigInt(n) >= k;
    }

    bool operator >= (BigInt n, BigInt k) {
        if(n.status == Status::Plus && k.status == Status::Minus) return true;
        if(n.status == Status::Minus && k.status == Status::Plus) return false;

        bool flag = true;
        if(n.status == Status::Minus) flag = false;

        if(n.digits > k.digits) return flag;
        if(n.digits < k.digits) return !flag;

        size_t i;
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
    
    template<typename T>
    bool operator < (T n, BigInt k) {
        return BigInt(n) < k;
    }

    bool operator < (BigInt n, BigInt k) {
        if(n.status == Status::Plus && k.status == Status::Minus) return false;
        if(n.status == Status::Minus && k.status == Status::Plus) return true;

        bool flag = true;
        if(n.status == Status::Minus) flag = false;

        if(n.digits < k.digits) return flag;
        if(n.digits > k.digits) return !flag;

        size_t i;
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

    template<typename T>
    bool operator <= (T n, BigInt k) {
        return BigInt(n) <= k;
    }

    bool operator <= (BigInt n, BigInt k) {
        if(n.status == Status::Plus && k.status == Status::Minus) return false;
        if(n.status == Status::Minus && k.status == Status::Plus) return true;

        bool flag = true;
        if(n.status == Status::Minus) flag = false;

        if(n.digits < k.digits) return flag;
        if(n.digits > k.digits) return !flag;

        size_t i;
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
    
    template<typename T>
    bool operator == (T n, BigInt k) {
        return BigInt(n) == k;
    }

    bool operator == (BigInt n, BigInt k) {
        if(n.status != k.status) return false;
        if(n.digits != k.digits) return false;
        for(size_t i = 0; i < n.digits; i++) {
            if(n.element[i] != k.element[i]) return false;
        }
        return true;
    }

    template<typename T>
    bool operator != (BigInt n, T k) {
        return n != BigInt(k);
    }

    template<typename T>
    bool operator != (T n, BigInt k) {
        return BigInt(n) != k;
    }

    bool operator != (BigInt n, BigInt k) {
        if(n.status != k.status) return true;
        if(n.digits != k.digits) return true;
        for(size_t i = 0; i < n.digits; i++) {
            if(n.element[i] != k.element[i]) return true;
        }
        return false;
    }

    // 論理演算
    BigInt operator ~ (BigInt n) {
        std::vector<size_t> b = n.to_binary();
        for(size_t i = 0; i < b.size(); i++) {
            b[i] ^= 1;
        }
        return binary_to_i(b);
    }

    template<typename T>
    BigInt operator & (BigInt n, T k) {
        return n & BigInt(k);
    }

    template<typename T>
    BigInt operator & (T n, BigInt k) {
        return BigInt(n) & k;
    }

    BigInt operator & (BigInt n, BigInt k) {
        if(n < k) return k & n;
        std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
        size_t i;
        for(i = 0; i < b_k.size(); i++) {
            b_k[i] &= b_n[i];
        }
        return binary_to_i(b_k);
    }

    template<typename T>
    BigInt operator | (BigInt n, T k) {
        return n | BigInt(k);
    }

    template<typename T>
    BigInt operator | (T n, BigInt k) {
        return BigInt(n) | k;
    }

    BigInt operator | (BigInt n, BigInt k) {
        if(n < k) return k | n;
        std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
        size_t i;
        for(i = 0; i < b_k.size(); i++) {
            b_n[i] |= b_k[i];
        }
        return binary_to_i(b_n);
    }

    template<typename T>
    BigInt operator ^ (BigInt n, T k) {
        return n ^ BigInt(k);
    }

    template<typename T>
    BigInt operator ^ (T n, BigInt k) {
        return BigInt(n) ^ k;
    }

    BigInt operator ^ (BigInt n, BigInt k) {
        if(n < k) return k ^ n;
        std::vector<size_t> b_n = n.to_binary(), b_k = k.to_binary();
        size_t i;
        for(i = 0; i < b_k.size(); i++) {
            b_n[i] ^= b_k[i];
        }
        return binary_to_i(b_n);
    }

    BigInt BigInt::operator << (size_t n) {
        if(n < 0) return *this >> (-n);
        std::vector<size_t> b = this->to_binary();
        for(size_t i = 0; i < n; ++i) {
            b.insert(b.begin(), 0);
        }
        return binary_to_i(b);
    }
    
    BigInt BigInt::operator >> (size_t n) {
        if(n < 0) return *this << (-n);
        std::vector<size_t> b = this->to_binary();
        auto begin = b.begin();
        if(n < b.size()) {
            b.erase(begin, begin + n);
            return binary_to_i(b);
        }
        return 0;
    }

    // I/O
    std::istream& operator >> (std::istream& in, BigInt& n) {
        std::string s;
        in >> s;
        n.digits = s.size();
        n.element.resize(n.digits + 50);
        for(size_t i = 0; i < n.digits; ++i) {
            if(s[n.digits - i - 1] == '-') n.status = Status::Minus;
            else n.element[i] = s[n.digits - i - 1] - '0';
        }
        return in;
    }

    std::ostream& operator << (std::ostream& out, const BigInt n) {
        if(n.status == Status::Minus) out << '-';
        for(int i = n.digits - 1; i >= 0; i--) {
            out << n.element[i];
        }
        return out;
    }
}

#endif
