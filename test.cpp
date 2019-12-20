#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include "BigInt.hpp"

void value_test(Rytar::BigInt n, Rytar::BigInt k) {
    std::cout << Rytar::BigInt("1234567890") << std::endl;
    std::cout << Rytar::BigInt(std::string("1234567890")) << std::endl;
    std::cout << Rytar::BigInt(1234567890) << std::endl;
    std::cout << n << " + " << k << " = " << (n + k) << std::endl;
    std::cout << n << " - " << k << " = " << (n - k) << std::endl;
    std::cout << n << " * " << k << " = " << (n * k) << std::endl;
    std::cout << n << " / " << k << " = " << (n / k) << std::endl;
    std::cout << n << " % " << k << " = " << (n % k) << std::endl;
    std::cout << n << " & " << k << " = " << (n & k) << std::endl;
    std::cout << n << " | " << k << " = " << (n | k) << std::endl;
    std::cout << n << " ^ " << k << " = " << (n ^ k) << std::endl;
    std::cout << n << " > " << k << " <=> " << (n > k) << std::endl;
    std::cout << n << " < " << k << " <=> " << (n < k) << std::endl;
    std::cout << n << " >> 2 = " << (n >> 2) << std::endl;
    std::cout << k << " << 3 = " << (k << 3) << std::endl;
    std::cout << (-n).abs() << " " << (-k).abs() << std::endl;
    std::cout << "char: " << n.to_char() << " " << k.to_char() << std::endl;
    std::cout << "unsigned char: " << n.to_unsigned_char() << " " << k.to_unsigned_char() << std::endl;
    std::cout << "short: " << n.to_short() << " " << k.to_short() << std::endl;
    std::cout << "unsigned short: " << n.to_unsigned_short() << " " << k.to_unsigned_short() << std::endl;
    std::cout << "int: " << n.to_int() << " " << k.to_int() << std::endl;
    std::cout << "unsigned int: " << n.to_unsigned_int() << " " << k.to_unsigned_int() << std::endl;
    std::cout << "long: " << n.to_long() << " " << k.to_long() << std::endl;
    std::cout << "unsigned long: " << n.to_unsigned_long() << " " << k.to_unsigned_long() << std::endl;
    std::cout << "long long: " << n.to_long_long() << " " << k.to_long_long() << std::endl;
    std::cout << "unsigned long long: " << n.to_unsigned_long_long() << " " << k.to_unsigned_long_long() << std::endl;
}

template<typename T>
void time_test(T a, T b) {
    T count = 1;
    auto start = std::chrono::system_clock::now();
    for(T i = a; i < b; i++) {
        // std::cout << i << " * " << i << " = " << (i * i) << std::endl;
        std::cout << i << " / " << count << " = " << (i / count) << std::endl;
        count++;
    }
    auto end = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>
              (end - start).count();
    std::cout << "time: " << ms << "[ms]" << std::endl;
    std::cout << "ave: " << (ms / count) << "[ms]" << std::endl;
}

template<typename T>
void time_test() {
    std::mt19937 mt;
    T a = mt(), b;
    auto start = std::chrono::system_clock::now();
    for(T i = 1; i < 100000; i++) {
        b = a + i;
        b = a - i;
        b = a * i;
        b = a / i;
    }
    auto end = std::chrono::system_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>
                (end - start).count();
    std::cout << a  << "^3 = " << (a * a * a) << std::endl;
    printf("ave: %.4f[ns]\n", (double)dur / 100000);
}

template<typename T>
void test() {
    T n = 16;
    auto start = std::chrono::system_clock::now();
    for(size_t i = 1; i < 4; i++) {
        std::cout << n << " * " << n << " = " << (n *= n) << std::endl;
        std::cout << n << " / " << i << " = " << (n / i) << std::endl;
    }
    auto end = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>
              (end - start).count();
    std::cout << "time: " << ms << "[μs]" << std::endl;
}


int main() {
    // Rytar::BigInt n, k;
    // std::cout << "input: ";
    // std::cin >> n >> k;
    time_test<Rytar::BigInt>();
    // value_test(n, k);

    // test<Rytar::BigInt>();
}
