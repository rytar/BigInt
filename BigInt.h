#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>
#include <string>

enum class Status {
    Plus = 1, Minus = -1
};

class BigInt {
private:
    std::vector<size_t> element;
    size_t digits;
    Status status;

public:
    BigInt();
    template<typename T>
    BigInt(T);
    BigInt(const std::string);

    // 算術演算
    BigInt operator + () const;
    BigInt operator - () const;
    BigInt& operator ++ ();
    BigInt operator ++ (int);
    BigInt& operator -- ();
    BigInt operator -- (int);

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
    friend BigInt karatsuba(BigInt, BigInt);
    template<typename T>
    friend BigInt operator * (BigInt, T);
    template<typename T>
    friend BigInt operator * (T, BigInt);
    friend BigInt operator * (BigInt, BigInt);
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

    // 関係演算子
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

    // 論理演算子
    BigInt& operator ~ ();
    template<typename T>
    friend BigInt operator & (BigInt, T);
    template<typename T>
    friend BigInt operator & (T, BigInt);
    friend BigInt operator & (BigInt, BigInt);
    friend BigInt operator | (BigInt, BigInt);
    friend BigInt operator ^ (BigInt, BigInt);
    template<typename T>
    BigInt operator << (T);
    template<typename T>
    BigInt operator >> (T);

    // 代入演算子

    BigInt& operator = (long long);
    BigInt& operator += (const BigInt);
    BigInt& operator -= (const BigInt);
    BigInt& operator *= (const BigInt);
    BigInt& operator /= (const BigInt);
    BigInt& operator %= (const BigInt);
    BigInt& operator &= (const BigInt);
    BigInt& operator |= (const BigInt);
    BigInt& operator ^= (const BigInt);
    BigInt& operator <<= (const size_t);
    BigInt& operator >>= (const size_t);

    // インデックス
    BigInt operator [] (const BigInt);

    // 入力・出力
    friend std::istream& operator >> (std::istream&, BigInt&);
    friend std::ostream& operator << (std::ostream&, const BigInt);

    // キャスト
    operator char() const;
    operator unsigned char() const;
    operator short() const;
    operator unsigned short() const;
    operator int() const;
    operator unsigned int() const;
    operator long() const;
    operator unsigned long() const;
    operator long long() const;
    operator unsigned long long() const;
    // operator float() const;
    // operator double() const;

    // その他
    size_t size();
    std::vector<size_t> to_binary();
    BigInt abs();
};

BigInt binary_to_i(std::vector<size_t>);

#endif
