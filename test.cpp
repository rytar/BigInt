#include <iostream>
#include <chrono>
#include <climits>
#include "BigInt.h"

void value_test(BigInt n, BigInt k) {
    std::cout << BigInt("1234567890") << std::endl;
    std::cout << BigInt(std::string("1234567890")) << std::endl;
    std::cout << BigInt(1234567890) << std::endl;
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
    std::cout << "char: " << (char)n << " " << (char)k << std::endl;
    std::cout << "unsigned char: " << (unsigned char)n << " " << (unsigned char)k << std::endl;
    std::cout << "short: " << (short)n << " " << (short)k << std::endl;
    std::cout << "unsigned short: " << (unsigned short)n << " " << (unsigned short)k << std::endl;
    std::cout << "int: " << (int)n << " " << (int)k << std::endl;
    std::cout << "unsigned int: " << (unsigned int)n << " " << (unsigned int)k << std::endl;
    std::cout << "long: " << (long)n << " " << (long)k << std::endl;
    std::cout << "unsigned long: " << (unsigned long)n << " " << (unsigned long)k << std::endl;
    std::cout << "long long: " << (long long)n << " " << (long long)k << std::endl;
    std::cout << "unsigned long long: " << (unsigned long long)n << " " << (unsigned long long)k << std::endl;
    // std::cout << "float: " << (float)n << " " << (float)k << std::endl;
    // std::cout << "double: " << (double)n << " " << (double)k << std::endl;
}

template<typename T>
void time_test(T a, T b) {
    T count = 1;
    auto start = std::chrono::system_clock::now();
    for(T i = a; i < b; i++) {
        std::cout << i << " * " << i << " = " << (i * i) << std::endl;
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
void test() {
    T n = 16;
    auto start = std::chrono::system_clock::now();
    for(size_t i = 1; i < 5; i++) {
        std::cout << n << " * " << n << " = " << (n *= n) << std::endl;
        std::cout << n << " / " << i << " = " << (n / i) << std::endl;
    }
    auto end = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>
              (end - start).count();
    std::cout << "time: " << ms << "[ms]" << std::endl;
}


int main() {
    BigInt n, k;
    // std::cout << "input: ";
    // std::cin >> n >> k;
    // time_test(n, k);
    // value_test(n, k);

    test<BigInt>();
}
