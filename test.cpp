#include <iostream>
#include <chrono>
#include <climits>
#include "BigInt.h"

template<typename T>
void value_test(T n, T k) {
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
}

template<typename T>
void time_test(T n, T k) {
    auto start = std::chrono::system_clock::now();
    for(T i = n; i < k; ++i)
        std::cout << i << " * " << i << " = " << (i * i) << std::endl;
    auto end = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>
              (end - start).count();
    std::cout << "time(ave): " << (ms / (k - n)) << "[ms]" << std::endl;
}

int main() {
    BigInt n, k;
    std::cout << "input: ";
    std::cin >> n >> k;
    time_test(n, k);
    value_test(n, k);
}
