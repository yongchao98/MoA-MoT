import math

def solve_last_digit():
    """
    This function calculates the last digit of A_k^B_k - B_k^A_k for a given k.
    The logic holds for any positive integer k. We use k=2 as an example.
    """
    k = 2  # Example value for k. The result is the same for any k >= 1.

    # A_k is a number of k+1 digits, all ones. (e.g., A_1=11, A_2=111)
    # B_k is 10^k.
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    print(f"Let's demonstrate with k = {k}:")
    print(f"A_k = (10^({k}+1) - 1) / 9 = {A_k}")
    print(f"B_k = 10^{k} = {B_k}")
    print("\nWe want to find the last digit of the expression A_k^B_k - B_k^A_k.")
    print("This is equivalent to (A_k^B_k - B_k^A_k) mod 10.")
    print("-" * 50)

    # Find the last digit of the first term, A_k^B_k
    # A number ending in 1, raised to any positive integer power, results in a number ending in 1.
    # We use pow(base, exp, mod) for efficient calculation of (base^exp) % mod.
    last_digit_term1 = pow(A_k, B_k, 10)
    print(f"The last digit of A_k^B_k (which is {A_k}^{B_k}) is {last_digit_term1}.")

    # Find the last digit of the second term, B_k^A_k
    # A number ending in 0 (for k>=1), raised to any positive integer power, results in a number ending in 0.
    # The exponent A_k is a positive integer.
    last_digit_term2 = pow(B_k, A_k, 10)
    print(f"The last digit of B_k^A_k (which is {B_k}^{A_k}) is {last_digit_term2}.")
    print("-" * 50)

    # The last digit of the difference is (last_digit_term1 - last_digit_term2) mod 10.
    final_result = (last_digit_term1 - last_digit_term2) % 10

    # As requested, printing each number in the final equation for the last digits.
    print("The final calculation for the last digit is based on the difference of the last digits:")
    print(f"{last_digit_term1} - {last_digit_term2} = {final_result}")
    print(f"\nTherefore, the last digit of A_k^B_k - B_k^A_k is {final_result}.")

solve_last_digit()