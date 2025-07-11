def solve_last_digit():
    """
    Calculates the last digit of A_k^(B_k) - B_k^(A_k) for any positive integer k.
    A_k = (10^(k+1) - 1) / 9
    B_k = 10^k
    """

    # Step 1: Find the last digit of A_k.
    # A_k is a number made of k+1 ones (e.g., 11, 111, ...).
    # Its last digit is always 1 for any positive k.
    last_digit_A_k = 1

    # Step 2: Find the last digit of B_k.
    # B_k is a power of 10 (e.g., 10, 100, ...).
    # Its last digit is always 0 for any positive k.
    last_digit_B_k = 0

    # Step 3: Find the last digit of A_k^B_k.
    # The last digit of a number ending in 1 raised to any positive integer power is always 1.
    # For k>=1, B_k is a positive integer.
    last_digit_term1 = 1

    # Step 4: Find the last digit of B_k^A_k.
    # The last digit of a number ending in 0 (except 0 itself) raised to any positive
    # integer power is always 0.
    # For k>=1, A_k is a positive integer.
    last_digit_term2 = 0

    # Step 5: Calculate the last digit of the difference.
    # The last digit of (Term1 - Term2) is (LastDigit(Term1) - LastDigit(Term2)) mod 10.
    final_last_digit = (last_digit_term1 - last_digit_term2) % 10

    # Print the explanation and the final equation with the computed last digits.
    print("This script determines the last digit of A_k^B_k - B_k^A_k.")
    print("The result is constant for any positive integer k.\n")
    print(f"The last digit of A_k is always {last_digit_A_k}.")
    print(f"The last digit of B_k is always {last_digit_B_k}.\n")
    print(f"The last digit of A_k^B_k, which is a number ending in {last_digit_A_k} raised to a positive integer power, is {last_digit_term1}.")
    print(f"The last digit of B_k^A_k, which is a number ending in {last_digit_B_k} raised to a positive integer power, is {last_digit_term2}.\n")

    print("The final equation for the last digit is:")
    print(f"(Last Digit of A_k^B_k) - (Last Digit of B_k^A_k) mod 10")
    print(f"= ({last_digit_term1} - {last_digit_term2}) mod 10")
    print(f"= {final_last_digit}")

# Run the function to display the result.
solve_last_digit()