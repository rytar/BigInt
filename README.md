# BigInt
This is a library for multiple‐precision arithmetic.
You can use overloaded operators dedicated to BigInt.

## Usage
Here is one sample code.
```
#include "BigInt.hpp"

using bint = Rytar::BigInt;

bint fib(bint n) {
    bint memo1 = 1, memo2 = 1, res = 0;
    for(bint i = 0; i < n; i++) {
        res = memo2;
        memo2 = memo1;
        memo1 = memo2 + res;
    }
    return res;
}

int main() {
    for(bint i = 0; i < 500; i++) std::cout << fib(i) << std::endl;
}
```
