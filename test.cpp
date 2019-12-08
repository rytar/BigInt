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
}

template<typename T>
void time_test(T n, T k) {
    auto start = std::chrono::system_clock::now();
    for(T i = n; i < k; ++i)
        std::cout << i << " * " << i << " = " << (i * i) << std::endl;
    auto end = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>
              (end - start).count();
    std::cout << "time: " << (ms / (k - n)) << "[ms]" << std::endl;
}

int main() {
    BigInt n, k;
    std::cout << "input: ";
    std::cin >> n >> k;
    time_test(n, k);
}