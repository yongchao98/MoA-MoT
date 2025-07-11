import math

def solve():
    """
    Calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given k.
    The problem states the last digit is constant for all positive integers k,
    so we can test with k=1.
    """
    k = 1

    # Definition of A_k and B_k
    # A_k = (10**(k+1) - 1) / 9 is a number with k+1 ones (e.g., 11, 111, ...)
    # B_k = 10**k
    Ak = (10**(k + 1) - 1) // 9
    Bk = 10**k

    # We need to find the last digit of Ak^Bk - Bk^Ak.
    # This is equivalent to (Ak^Bk - Bk^Ak) % 10.

    # Find the last digit of the first term, Ak^Bk
    # The last digit of a power can be found using modular exponentiation.
    # pow(base, exp, mod) is equivalent to (base ** exp) % mod
    # Note: For Ak, the last digit is always 1. Any power of a number ending in 1 also ends in 1.
    last_digit_A_pow_B = pow(Ak, Bk, 10)

    # Find the last digit of the second term, Bk^Ak
    # Note: For Bk (k>=1), the last digit is always 0. Any positive power of a number ending in 0
    # also ends in 0.
    last_digit_B_pow_A = pow(Bk, Ak, 10)

    # Calculate the last digit of the difference.
    # The ( ... + 10) % 10 is a standard way to handle potential negative
    # results in modular arithmetic, ensuring the result is always positive.
    final_last_digit = (last_digit_A_pow_B - last_digit_B_pow_A + 10) % 10

    # Print the explanation and the final equation with numbers
    print(f"Let's test for k = {k}.")
    print(f"A_k = (10^({k}+1) - 1) / 9  =>  A_{k} = {Ak}")
    print(f"B_k = 10^{k}  =>  B_{k} = {Bk}")
    print(f"We need the last digit of the expression: {Ak}^{Bk} - {Bk}^{Ak}")
    print(f"Last digit of the first term ({Ak}^{Bk}): {last_digit_A_pow_B}")
    print(f"Last digit of the second term ({Bk}^{Ak}): {last_digit_B_pow_A}")
    print(f"The last digit of the difference is ({last_digit_A_pow_B} - {last_digit_B_pow_A}) mod 10, which is {final_last_digit}.")

solve()
<<<1>>>