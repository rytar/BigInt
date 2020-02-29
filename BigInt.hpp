#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>
#include <unordered_map>
#include <utility>

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
            element[0] = 0;
            digits = 1;
            status = Status::Plus;
        }

        template<typename T>
        BigInt(T n) {
            if(n == 0) {
                element.resize(50);
                element[0] = 0;
                digits = 1;
                status = Status::Plus;
                return;
            }

            if(n < 0) {
                status = Status::Minus;
                n = -n;
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
            for(size_t i = 0; i < digits - 1; ++i) {
                element[i] = s[digits - i - 1] - '0';
            }
            if(s[0] == '-') {
                digits--;
                status = Status::Minus;
            }
            else element[digits - 1] = s[0] - '0';
        }
        
        BigInt(const std::string s) {
            digits = s.size();
            element.resize(digits + 50);
            status = Status::Plus;
            for(size_t i = 0; i < digits - 1; ++i) {
                element[i] = s[digits - i - 1] - '0';
            }
            if(s[0] == '-') {
                digits--;
                status = Status::Minus;
            }
            else element[digits - 1] = s[0] - '0';
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
            BigInt n1 = 0, n2 = 0, k1 = 0, k2 = 0, r2 = 0, r1 = 0, r0 = 0;
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
            r2 = n2 * k2;
            r0 = n1 * k1;
            return r0 + (r2 + r0 - (n2 - n1) * (k2 - k1)) * B + r2 * B * B;
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
            if(*this > CHAR_MAX || *this < CHAR_MIN + 1) return 0;
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
            if(*this > SHRT_MAX || *this < SHRT_MIN + 1) return 0;
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
            if(*this > INT_MAX || *this < INT_MIN + 1) return 0;
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
            if(*this > LONG_MAX || *this < LONG_MIN + 1) return 0;
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
            if(*this > LLONG_MAX || *this < LLONG_MIN + 1) return 0;
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

        std::string to_string() const {
            std::string s(this->digits, 'a');
            for(size_t i = 0; i < this->digits; i++) {
                s[this->digits - i - 1] = this->element[i] + '0';
            }
            return s;
        }

        // その他
        size_t size() {
            return digits;
        }

        BigInt abs() {
            BigInt n = *this;
            if(n.status == Status::Minus) n.status = Status::Plus;
            return n;
        }

        std::vector<size_t> to_binary() {
            std::vector<size_t> binary;
            if(*this == 0) {
                binary.push_back(0);
                return binary;
            }
            BigInt n = *this;
            while(n != 0) {
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

        static std::unordered_map<std::string, BigInt> add;
        std::string p = n.to_string() + "+" + k.to_string();
        try {
            return add.at(p);
        }
        catch(std::out_of_range) {
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
            return (add[p] = n);
        }
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
        static std::unordered_map<std::string, BigInt> sub;
        std::string p = n.to_string() + "-" + k.to_string();
        try {
            return sub.at(p);
        }
        catch(std::out_of_range) {
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
            return (sub[p] = n);
        }
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
        if(n == 1) return k;
        if(k == 1) return n;

        static std::unordered_map<std::string, BigInt> mul;
        std::string p = n.to_string() + "*" + k.to_string();
        try {
            return mul.at(p);
        }
        catch(std::out_of_range) {
            size_t digits_sum = n.digits + k.digits;
            auto max = std::pow(10, digits_sum);
            if(max < CHAR_MAX) return (mul[p] = n.to_char() * k.to_char());
            if(max < SHRT_MAX) return (mul[p] = n.to_short() * k.to_short());
            if(max < INT_MAX) return (mul[p] = n.to_int() * k.to_int());
            if(max < LONG_MAX) return (mul[p] = n.to_long() * k.to_long());
            if(max < LLONG_MAX) return (mul[p] = n.to_long_long() * k.to_long_long());

            size_t i, j;
            BigInt tmp = k;

            if(n.status == k.status)
                k.status = Status::Plus;
            else k.status = Status::Minus;

            k.element.resize(k.element.size() + 10);

            if(n.digits / 2 < k.digits && n.digits > 4) {
                size_t count;
                for(count = 0; k.element[count] == 0; count++);
                std::vector<size_t> add(count, 0);
                n.element.insert(n.element.begin(), add.begin(), add.end());
                k.element.erase(k.element.begin(), k.element.begin() + count);
                n.digits += count;
                k.digits -= count;
                
                if(k.digits > 2) {
                    return (mul[p] = BigInt::karatsuba(n, k));
                }
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
            return (mul[p] = k);
        }
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
        if(n < 0) return -(-n / k);
        if(k < 0) return -(n / -k);
        if(n < k) return 0;
        if(n == k) return 1;
        static std::unordered_map<std::string, BigInt> div;
        std::string p = n.to_string() + "/" + k.to_string();
        try {
            return div[p];
        }
        catch(std::out_of_range) {
            if(n < CHAR_MAX) return (div[p] = n.to_char() / k.to_char());
            if(n < SHRT_MAX) return (div[p] = n.to_short() / k.to_short());
            if(n < INT_MAX)  return (div[p] = n.to_int() / k.to_int());
            if(n < LONG_MAX) return (div[p] = n.to_long() / k.to_long());
            if(n < LLONG_MAX) return (div[p] = n.to_long_long() / k.to_long_long());
            if(n.digits > 2 * k.digits + 2) {
                return (div[p] = BigInt::karatsuba_for_div(n, k));
            }
            BigInt q = 0, add = 0, next = 1;
            while(n - q < k) {
                if(k * (q + next) <= 0) {
                    add = next;
                    next *= k;
                }
                else {
                    q += add;
                    add = 0;
                    next = 1;
                }
            }
            return (div[p] = q);
        }
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
        if(n < k) return n;
        static std::unordered_map<std::string, BigInt> mod;
        std::string p = n.to_string() + "%" + k.to_string();
        try {
            return mod.at(p);
        }
        catch(std::out_of_range) {
            return (mod[p] = n - k * (n / k));
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status == Status::Plus && k.status == Status::Minus) return true;
        if(n.status == Status::Minus && k.status == Status::Plus) return false;
        if(n == k) return false;

        static std::unordered_map<std::string, bool> gt;
        std::string p = n.to_string() + ">" + k.to_string();
        try {
            return gt.at(p);
        }
        catch(std::out_of_range) {
            gt[p] = true;
            if(n.status == Status::Minus) gt[p] = false;

            if(n.digits > k.digits) return gt[p];
            if(n.digits < k.digits) return (gt[p] = !gt[p]);

            size_t i;
            for(i = n.digits - 1; i >= 0; i--) {
                if(n.element[i] > k.element[i]) return gt[p];
                if(n.element[i] < k.element[i]) return (gt[p] = !gt[p]);
            }
            return (gt[p] = false);
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status == Status::Plus && k.status == Status::Minus) return true;
        if(n.status == Status::Minus && k.status == Status::Plus) return false;
        if(n == k) return true;

        static std::unordered_map<std::string, bool> gte;
        std::string p = n.to_string() + ">=" + k.to_string();
        try {
            return gte.at(p);
        }
        catch(std::out_of_range) {
            gte[p] = true;
            if(n.status == Status::Minus) gte[p] = false;

            if(n.digits > k.digits) return gte[p];
            if(n.digits < k.digits) return (gte[p] = !gte[p]);

            size_t i;
            for(i = n.digits - 1; i >= 0; i--) {
                if(n.element[i] > k.element[i]) return gte[p];
                if(n.element[i] < k.element[i]) return (gte[p] = !gte[p]);
            }
            return (gte[p] = true);
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status == Status::Plus && k.status == Status::Minus) return false;
        if(n.status == Status::Minus && k.status == Status::Plus) return true;
        if(n == k) return false;

        static std::unordered_map<std::string, bool> lt;
        std::string p = n.to_string() + "<" + k.to_string();
        try {
            return lt.at(p);
        }
        catch(std::out_of_range) {
            lt[p] = true;
            if(n.status == Status::Minus) lt[p] = false;

            if(n.digits < k.digits) return lt[p];
            if(n.digits > k.digits) return (lt[p] = !lt[p]);

            size_t i;
            for(i = n.digits - 1; i >= 0; i--) {
                if(n.element[i] < k.element[i]) return lt[p];
                if(n.element[i] > k.element[i]) return (lt[p] = !lt[p]);
            }
            return (lt[p] = false);
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status == Status::Plus && k.status == Status::Minus) return false;
        if(n.status == Status::Minus && k.status == Status::Plus) return true;
        if(n == k) return true;

        static std::unordered_map<std::string, bool> lte;
        std::string p = n.to_string() + "<=" + k.to_string();
        try {
            return lte.at(p);
        }
        catch(std::out_of_range) {
            lte[p] = true;
            if(n.status == Status::Minus) lte[p] = false;

            if(n.digits < k.digits) return lte[p];
            if(n.digits > k.digits) return (lte[p] = !lte[p]);

            size_t i;
            for(i = n.digits - 1; i >= 0; i--) {
                if(n.element[i] < k.element[i]) return lte[p];
                if(n.element[i] > k.element[i]) return (lte[p] = !lte[p]);
            }
            return (lte[p] = true);
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status != k.status) return false;
        if(n.digits != k.digits) return false;
        static std::unordered_map<std::string, bool> eql;
        std::string p = n.to_string() + "==" + k.to_string();
        try {
            return eql.at(p);
        }
        catch(std::out_of_range) {
            for(size_t i = 0; i < n.digits; i++) {
                if(n.element[i] != k.element[i]) return (eql[p] = false);
            }
            return (eql[p] = true);
        }
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
        if(n.digits == 0 || k.digits == 0) std::cerr << "digits 0 error" << std::endl;
        if(n.status != k.status) return true;
        if(n.digits != k.digits) return true;
        static std::unordered_map<std::string, bool> neq;
        std::string p = n.to_string() + "==" + k.to_string();
        try {
            return neq.at(p);
        }
        catch(std::out_of_range) {
            for(size_t i = 0; i < n.digits; i++) {
                if(n.element[i] != k.element[i]) return (neq[p] = true);
            }
            return (neq[p] = false);
        }
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
        BigInt abs = this->abs();
        std::vector<size_t> b = abs.to_binary();
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
            if(s[n.digits - i - 1] == '-') {
                n.digits--;
                n.status = Status::Minus;
            }
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
