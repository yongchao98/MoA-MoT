def find_last_digit():
    """
    Calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given positive integer k.
    The problem states the last digit is the same for all k, so we test for an arbitrary k, k=1.
    """
    k = 1

    # Define A_k and B_k.
    # A_k = (10^(k+1) - 1) / 9 is a repunit with k+1 digits.
    # B_k = 10^k.
    # We don't need the full values, only their last digits, for the calculation.
    # However, we need the full values for the exponents in the pow function.
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    # We need to find the last digit of the expression, which is the expression modulo 10.
    # We can use modular exponentiation to find the last digit of each term.
    # Python's built-in pow(base, exp, mod) is efficient for this purpose.

    # Last digit of the first term: A_k^B_k mod 10
    last_digit_term1 = pow(A_k, B_k, 10)

    # Last digit of the second term: B_k^A_k mod 10
    last_digit_term2 = pow(B_k, A_k, 10)

    # The last digit of the difference (A_k^B_k - B_k^A_k).
    # The final result is (last_digit_term1 - last_digit_term2) mod 10.
    # Python's % operator handles negative results correctly for finding the last digit.
    # For example, (3 - 7) % 10 = -4 % 10 = 6.
    final_last_digit = (last_digit_term1 - last_digit_term2) % 10

    # Output the numbers in the final equation for the last digit, as requested.
    print(f"The last digit of the first term, A_k^B_k, is {last_digit_term1}.")
    print(f"The last digit of the second term, B_k^A_k, is {last_digit_term2}.")
    print("The final equation for the last digit of the expression is:")
    print(f"({last_digit_term1} - {last_digit_term2}) mod 10 = {final_last_digit}")
    print(f"\nThe constant last digit is {final_last_digit}.")

find_last_digit()