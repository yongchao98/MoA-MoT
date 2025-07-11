import math

def solve_last_digit():
    """
    This function calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given k.
    The problem states the last digit is the same for all positive integers k.
    We will use k=2 as an example.
    """
    k = 2

    # Define A_k and B_k
    # A_k = (10^(k+1) - 1) / 9 is a number with k+1 ones (e.g., 11, 111, ...).
    # B_k = 10^k (e.g., 10, 100, ...).
    Ak = (10**(k + 1) - 1) // 9
    Bk = 10**k

    print(f"Let's demonstrate the calculation for k = {k}:")
    print(f"A_k = {Ak}")
    print(f"B_k = {Bk}")
    print(f"We want to find the last digit of the expression: {Ak}^{Bk} - {Bk}^{Ak}")
    print("-" * 30)

    # To find the last digit, we compute the value modulo 10.
    # First term: A_k^(B_k)
    # The last digit of Ak is always 1. Any power of 1 is 1.
    last_digit_Ak_pow_Bk = pow(Ak, Bk, 10)
    print(f"The last digit of the first term, {Ak}^{Bk}, is {last_digit_Ak_pow_Bk}.")

    # Second term: B_k^(A_k)
    # The last digit of Bk is always 0 (for k>=1). Any positive power of 0 is 0.
    # Ak is always a positive integer for k>=1.
    last_digit_Bk_pow_Ak = pow(Bk, Ak, 10)
    print(f"The last digit of the second term, {Bk}^{Ak}, is {last_digit_Bk_pow_Ak}.")
    print("-" * 30)

    # The last digit of the difference is (last_digit_1 - last_digit_2) mod 10.
    # We add 10 before the modulo to handle potential negative results correctly.
    final_last_digit = (last_digit_Ak_pow_Bk - last_digit_Bk_pow_Ak + 10) % 10

    # Output the final equation for the last digits
    print("The final calculation for the last digit is:")
    print(f"({last_digit_Ak_pow_Bk} - {last_digit_Bk_pow_Ak}) mod 10 = {final_last_digit}")
    
    print("\nThis holds true for any positive integer k.")
    print(f"The constant last digit is: {final_last_digit}")


solve_last_digit()