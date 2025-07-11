import math

def solve_last_digit():
    """
    Calculates the last digit of A_k^{B_k} - B_k^{A_k}.

    A_k = (10^(k+1) - 1) / 9, which is a sequence of k+1 ones.
    B_k = 10^k.

    The last digit of a number N is N % 10.
    The last digit of the expression is (A_k^{B_k} - B_k^{A_k}) % 10.
    This can be calculated as ((A_k % 10)^(B_k) - (B_k % 10)^(A_k)) % 10,
    with careful application of modular arithmetic rules.
    """
    # The logic is independent of k for any positive integer k.
    # We use k=2 as a concrete example to demonstrate the calculation.
    k = 2

    print(f"Let's solve for k = {k}. The logic holds for any positive integer k.")
    print("-" * 50)

    # 1. Determine the last digit of A_k
    # A_k is a number made of k+1 ones (e.g., 11, 111, ...).
    # So, its last digit is always 1.
    last_digit_Ak = 1
    a_k_example = (10**(k + 1) - 1) // 9
    print(f"For k={k}, A_k = {a_k_example}.")
    print(f"The last digit of A_k is always {last_digit_Ak}.")
    print("-" * 50)

    # 2. Determine the last digit of B_k
    # B_k = 10^k is always a multiple of 10 for k >= 1.
    # So, its last digit is always 0.
    last_digit_Bk = 0
    b_k_example = 10**k
    print(f"For k={k}, B_k = {b_k_example}.")
    print(f"The last digit of B_k is always {last_digit_Bk}.")
    print("-" * 50)

    # 3. Determine the last digit of the first term: A_k^{B_k}
    # The last digit of A_k is 1. Any positive integer power of a number
    # ending in 1 also ends in 1. The exponent B_k is a positive integer.
    # We can use pow(base, exp, mod) for modular exponentiation.
    # Note: the exponent B_k can be a large number, but we only need its effect.
    last_digit_term1 = 1
    print(f"The last digit of A_k^(B_k) is the last digit of {last_digit_Ak}^({b_k_example}).")
    print(f"Since the base's last digit is 1, the result's last digit is {last_digit_term1}.")
    print("-" * 50)

    # 4. Determine the last digit of the second term: B_k^{A_k}
    # The last digit of B_k is 0. Any positive integer power of a number
    # ending in 0 also ends in 0. The exponent A_k is a positive integer.
    last_digit_term2 = 0
    print(f"The last digit of B_k^(A_k) is the last digit of {last_digit_Bk}^({a_k_example}).")
    print(f"Since the base's last digit is 0 and the exponent is positive, the result's last digit is {last_digit_term2}.")
    print("-" * 50)

    # 5. Calculate the last digit of the difference
    # result = (last_digit_term1 - last_digit_term2) mod 10
    # Add 10 before the modulo to handle potential negative results correctly.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    print("The final calculation for the last digit is:")
    print(f"(Last digit of A_k^(B_k) - Last digit of B_k^(A_k)) mod 10")
    print("\nHere is the final equation with the numbers plugged in:")
    print(f"({last_digit_term1} - {last_digit_term2} + 10) % 10 = {final_last_digit}")

solve_last_digit()