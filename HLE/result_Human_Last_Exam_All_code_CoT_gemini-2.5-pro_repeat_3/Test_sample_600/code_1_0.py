def solve_last_digit():
    """
    This function finds the last digit of A_k^B_k - B_k^A_k for any positive integer k.
    The logic is based on modular arithmetic (finding the value modulo 10).
    The script demonstrates the calculation for a specific k to illustrate the process.
    """
    # We can choose any positive integer for k, the result will be the same. Let's use k=2.
    k = 2

    print(f"Let's solve the problem for a sample value of k = {k}\n")

    # Definitions:
    # A_k = (10**(k+1) - 1) / 9 is a number with k+1 digits, all of them 1s.
    # B_k = 10**k

    # For large k, the actual numbers are enormous. We only need their last digits for the calculation.

    print("Step 1: Determine the last digit of A_k and B_k.")
    
    # For any k >= 1, A_k is a number like 11, 111, 1111, ...
    # So, its last digit is always 1.
    a_k_last_digit = 1
    
    # For any k >= 1, B_k is a number like 10, 100, 1000, ...
    # So, its last digit is always 0.
    b_k_last_digit = 0
    
    # We'll calculate the actual values of A_k and B_k for k=2 just for demonstration.
    A_k = (10**(k+1) - 1) // 9
    B_k = 10**k
    
    print(f"For k={k}:")
    print(f"A_k = (10^({k}+1) - 1) / 9 = {A_k}")
    print(f"B_k = 10^{k} = {B_k}")
    print(f"The last digit of A_{k} is {a_k_last_digit}.")
    print(f"The last digit of B_{k} is {b_k_last_digit}.\n")

    print("Step 2: Find the last digit of the first term, A_k^B_k.")
    # The last digit of a power depends on the last digit of the base.
    # The last digit of A_k is 1. Any positive integer power of a number ending in 1 also ends in 1.
    # We use pow(base, exp, mod) for modular exponentiation.
    term1_last_digit = pow(a_k_last_digit, B_k, 10)
    print(f"The last digit of A_{k}^{{B_k}} is equivalent to {a_k_last_digit}^{{{B_k}}} mod 10, which is {term1_last_digit}.\n")

    print("Step 3: Find the last digit of the second term, B_k^A_k.")
    # The last digit of B_k is 0. Any positive integer power of a number ending in 0 also ends in 0.
    # A_k is a positive integer for k>=1.
    term2_last_digit = pow(b_k_last_digit, A_k, 10)
    print(f"The last digit of B_{k}^{{A_k}} is equivalent to {b_k_last_digit}^{{{A_k}}} mod 10, which is {term2_last_digit}.\n")

    print("Step 4: Find the last digit of the final expression A_k^B_k - B_k^A_k.")
    # The last digit of a difference (X - Y) is (last_digit(X) - last_digit(Y)) mod 10.
    # We add 10 before taking the modulo to ensure the result is non-negative.
    final_last_digit = (term1_last_digit - term2_last_digit + 10) % 10
    
    print("The final equation for the last digit is:")
    print(f"LastDigit(A_k^B_k - B_k^A_k) = (LastDigit(A_k^B_k) - LastDigit(B_k^A_k)) mod 10")
    print(f"Result = ({term1_last_digit} - {term2_last_digit} + 10) % 10 = {final_last_digit}")
    
    print(f"\nThus, for any positive integer k, the last digit of the expression is always {final_last_digit}.")

solve_last_digit()