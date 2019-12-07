#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>

enum class Status {
    Plus, Minus
};

class BigInt {
private:
    std::vector<size_t> element;
    size_t digits;
    Status status;

public:
    BigInt();
    BigInt(long long);

    // 算術演算
    BigInt operator + () const;
    BigInt operator - () const;
    BigInt& operator ++ ();
    BigInt operator ++ (int);
    BigInt& operator -- ();
    BigInt operator -- (int);

    friend BigInt operator + (BigInt, BigInt);
    friend BigInt operator - (BigInt, BigInt);
    friend BigInt karatsuba(BigInt, BigInt);
    friend BigInt operator * (BigInt, BigInt);
    friend BigInt operator / (BigInt, BigInt);
    friend BigInt operator % (BigInt, BigInt);

    // 関係演算子
    friend bool operator > (BigInt, BigInt);
    friend bool operator >= (BigInt, BigInt);
    friend bool operator < (BigInt, BigInt);
    friend bool operator <= (BigInt, BigInt);
    friend bool operator == (BigInt, BigInt);
    friend bool operator != (BigInt, BigInt);

    // 論理演算子
    BigInt& operator ~ ();
    friend BigInt operator & (BigInt, BigInt);
    friend BigInt operator | (BigInt, BigInt);
    friend BigInt operator ^ (BigInt, BigInt);
    BigInt& operator << (BigInt);
    BigInt& operator >> (BigInt);

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
    BigInt& operator <<= (const BigInt);
    BigInt& operator >>= (const BigInt);

    // インデックス
    BigInt operator [] (const BigInt);

    // 入力・出力
    friend std::istream& operator >> (std::istream&, BigInt&);
    friend std::ostream& operator << (std::ostream&, const BigInt);

    // その他
    size_t size();
    std::vector<size_t> to_binary();
};

BigInt binary_to_i(std::vector<size_t>);

#endif