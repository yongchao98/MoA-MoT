def solve_last_digit():
    """
    This function determines the last digit of A_k^(B_k) - B_k^(A_k)
    for any positive integer k.
    """

    print("The problem asks for the last digit of the expression A_k^(B_k) - B_k^(A_k).")
    print("This last digit is constant for any positive integer k.")
    print("We can find this by analyzing the last digits of each part of the expression.\n")

    # Step 1: Find the last digit of A_k
    # A_k = (10**(k+1) - 1) / 9 results in a number made of k+1 ones (a repunit).
    # For example, A_1 = 11, A_2 = 111.
    last_digit_A = 1
    print(f"1. For any positive integer k, A_k is a number like 11, 111, ..., so its last digit is always {last_digit_A}.")

    # Step 2: Find the last digit of B_k
    # B_k = 10**k. For any positive integer k, B_k is a multiple of 10.
    # For example, B_1 = 10, B_2 = 100.
    last_digit_B = 0
    print(f"2. For any positive integer k, B_k is a power of 10, so its last digit is always {last_digit_B}.")

    # Step 3: Find the last digit of A_k^(B_k)
    # The base ends in 1. Any positive integer power of a number ending in 1 also ends in 1.
    last_digit_A_pow_B = 1
    print(f"3. The last digit of A_k^(B_k) is based on (...{last_digit_A})^B_k, which is {last_digit_A_pow_B}.")

    # Step 4: Find the last digit of B_k^(A_k)
    # The base ends in 0. Any positive integer power of a number ending in 0 also ends in 0.
    last_digit_B_pow_A = 0
    print(f"4. The last digit of B_k^(A_k) is based on (...{last_digit_B})^A_k, which is {last_digit_B_pow_A}.")

    # Step 5: Calculate the final last digit
    final_last_digit = (last_digit_A_pow_B - last_digit_B_pow_A + 10) % 10
    print("\nFinally, we find the last digit of the difference.")
    print("The equation using the last digits of each term is:")
    print(f"{last_digit_A_pow_B} - {last_digit_B_pow_A} = {final_last_digit}")
    
    print(f"\nThus, the last digit of the expression for any positive integer k is {final_last_digit}.")


solve_last_digit()