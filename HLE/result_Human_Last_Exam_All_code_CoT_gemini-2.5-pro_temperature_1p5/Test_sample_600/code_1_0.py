import math

def solve_problem():
    """
    This function calculates the last digit of A_k^B_k - B_k^A_k
    for a given k, and prints the step-by-step reasoning.
    """
    # Let's use k=1 for demonstration, as the last digit is the same for all k > 0.
    k = 1

    # Calculate A_k and B_k
    # A_k = (10**(k+1) - 1) / 9 is a repunit number (all 1s)
    # B_k = 10**k
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    print(f"Let's analyze the expression for k = {k}.")
    print(f"A_k = (10^({k}+1) - 1) / 9  =>  A_{k} = {A_k}")
    print(f"B_k = 10^{k}  =>  B_{k} = {B_k}")
    print(f"We want to find the last digit of A_k^B_k - B_k^A_k, which is ({A_k}^{B_k} - {B_k}^{A_k}) mod 10.\n")

    # Step 1: Find the last digit of A_k^B_k
    last_digit_A = A_k % 10
    # The last digit of a number ending in 1, raised to any integer power, is 1.
    last_digit_A_pow_B = pow(last_digit_A, B_k, 10)
    print(f"First, we find the last digit of {A_k}^{B_k}.")
    print(f"The last digit of A_{k} ({A_k}) is {last_digit_A}.")
    print(f"Any integer power of a number ending in 1 also ends in 1.")
    print(f"So, the last digit of {A_k}^{B_k} is {last_digit_A_pow_B}.\n")

    # Step 2: Find the last digit of B_k^A_k
    last_digit_B = B_k % 10
    # The last digit of a number ending in 0, raised to any positive integer power, is 0.
    # Since A_k is always a positive integer, this holds.
    last_digit_B_pow_A = 0
    if B_k > 0: # Check to ensure B_k is not 0, though for k>=1 it is always a power of 10.
        last_digit_B_pow_A = pow(last_digit_B, A_k, 10) if A_k > 0 else 1 # 0^0 is 1 by convention here

    print(f"Next, we find the last digit of {B_k}^{A_k}.")
    print(f"The last digit of B_{k} ({B_k}) is {last_digit_B}.")
    print(f"Any positive integer power of a number ending in 0 also ends in 0.")
    print(f"So, the last digit of {B_k}^{A_k} is {last_digit_B_pow_A}.\n")

    # Step 3: Calculate the final last digit
    # Use (a - b + 10) % 10 to correctly handle modular arithmetic with subtraction
    final_last_digit = (last_digit_A_pow_B - last_digit_B_pow_A + 10) % 10
    print("Finally, we calculate the last digit of the difference:")
    print(f"(last digit of {A_k}^{B_k}) - (last digit of {B_k}^{A_k}) mod 10")
    print(f"= ({last_digit_A_pow_B} - {last_digit_B_pow_A}) mod 10")

    print("\nThe final equation is:")
    print(f"{last_digit_A_pow_B} - {last_digit_B_pow_A} = {final_last_digit}")


solve_problem()