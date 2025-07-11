import math

def solve_last_digit():
    """
    This function calculates the last digit of A_k^B_k - B_k^A_k for a given k.
    A_k = (10**(k+1) - 1) / 9
    B_k = 10**k
    """
    # We can test this for any positive integer k. The result should be the same.
    k = 3

    print(f"The problem is to find the last digit of A_k^B_k - B_k^A_k for any positive integer k.")
    print(f"Let's demonstrate the solution for k = {k}.")
    print("-" * 30)

    # Calculate A_k and B_k
    # A_k is a sequence of k+1 ones (e.g., k=1 -> 11, k=2 -> 111)
    A_k = (10**(k + 1) - 1) // 9
    # B_k is 10 to the power of k
    B_k = 10**k

    print(f"Step 1: Calculate A_k and B_k for k={k}")
    print(f"A_{k} = (10^({k}+1) - 1) / 9 = {A_k}")
    print(f"B_{k} = 10^{k} = {B_k}")
    print()

    print(f"The expression we need to evaluate the last digit of is: {A_k}^{B_k} - {B_k}^{A_k}")
    print()
    
    # Calculate the last digit of the first term: A_k^B_k
    # We use modular exponentiation pow(base, exp, mod) which is efficient.
    # The last digit of A_k is A_k % 10, which is always 1.
    # Any power of a number ending in 1 results in a number ending in 1.
    last_digit_term1 = pow(A_k, B_k, 10)
    
    print(f"Step 2: Find the last digit of the first term, {A_k}^{B_k}")
    print(f"The last digit is pow({A_k}, {B_k}, 10) = {last_digit_term1}")
    print()

    # Calculate the last digit of the second term: B_k^A_k
    # The last digit of B_k is B_k % 10, which is always 0 for k>=1.
    # Any positive power of a number ending in 0 results in a number ending in 0.
    # A_k is always a positive integer, so the exponent is > 0.
    last_digit_term2 = pow(B_k, A_k, 10)

    print(f"Step 3: Find the last digit of the second term, {B_k}^{A_k}")
    print(f"The last digit is pow({B_k}, {A_k}, 10) = {last_digit_term2}")
    print()

    # The final last digit is (last_digit_term1 - last_digit_term2) mod 10
    # The formula (a - b + m) % m handles potential negative results from subtraction.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    print(f"Step 4: Calculate the last digit of the final expression")
    print(f"The final calculation is ({last_digit_term1} - {last_digit_term2}) mod 10")
    print(f"The constant last digit is: {final_last_digit}")

solve_last_digit()