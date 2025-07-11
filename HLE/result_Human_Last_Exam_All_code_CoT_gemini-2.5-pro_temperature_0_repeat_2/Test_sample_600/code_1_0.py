def find_last_digit():
    """
    This function calculates the last digit of the expression A_k^(B_k) - B_k^(A_k).
    The problem states the last digit is constant for any positive integer k.
    We will demonstrate this by choosing k=1.
    """
    k = 1

    # For any positive integer k:
    # A_k = (10**(k+1) - 1) / 9 is a number like 11, 111, 1111, ...
    # B_k = 10**k is a number like 10, 100, 1000, ...

    # Let's calculate A_k and B_k for our chosen k.
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    print(f"The problem is to find the last digit of A_k^B_k - B_k^A_k.")
    print(f"Let's test for k = {k}.")
    print(f"A_{k} = {A_k}")
    print(f"B_{k} = {B_k}")
    print(f"The expression is: {A_k}^{B_k} - {B_k}^{A_k}")
    print("-" * 20)

    # To find the last digit of a number, we calculate the number modulo 10.
    # We can use the property that (X - Y) mod N = ((X mod N) - (Y mod N)) mod N.

    # Find the last digit of the first term, A_k^B_k.
    # The last digit of A_k is always 1. Any power of a number ending in 1 also ends in 1.
    last_digit_term1 = pow(A_k, B_k, 10)

    # Find the last digit of the second term, B_k^A_k.
    # The last digit of B_k (for k>=1) is 0. Any positive power of a number ending in 0 also ends in 0.
    last_digit_term2 = pow(B_k, A_k, 10)

    # The final equation for the last digits is the difference of the last digits.
    # We add 10 before the final modulo to ensure the result is non-negative.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    print(f"The last digit of {A_k}^{B_k} is {last_digit_term1}.")
    print(f"The last digit of {B_k}^{A_k} is {last_digit_term2}.")
    print("\nThe final equation for the last digit is the last digit of the difference:")
    print(f"Last digit of ({last_digit_term1} - {last_digit_term2})")
    print(f"The result is {final_last_digit}.")

find_last_digit()