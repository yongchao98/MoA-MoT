import math

def solve_last_digit():
    """
    Calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given k.
    The problem states the last digit is the same for every positive integer k.
    We'll demonstrate with k=1.
    """
    k = 1

    # Define A_k and B_k
    # A_k = (10^(k+1) - 1) / 9
    # B_k = 10^k
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    # We need to find the last digit of A_k^B_k - B_k^A_k.
    # This is equivalent to (A_k^B_k - B_k^A_k) mod 10.
    # We can use the pow(base, exp, mod) function for efficiency.

    # Find the last digit of the first term: A_k^B_k
    # Since A_k's last digit is 1, any power of A_k will have a last digit of 1.
    last_digit_A_k_term = pow(A_k, B_k, 10)

    # Find the last digit of the second term: B_k^A_k
    # Since B_k is a multiple of 10, any positive power of B_k will have a last digit of 0.
    last_digit_B_k_term = pow(B_k, A_k, 10)

    # The last digit of the difference.
    # (a - b) % 10. Python's % operator handles negative results correctly for this purpose.
    final_last_digit = (last_digit_A_k_term - last_digit_B_k_term) % 10

    # Print the breakdown of the calculation and the final equation.
    print(f"Let's test for k = {k}:")
    print(f"A_{k} = {A_k}")
    print(f"B_{k} = {B_k}")
    print(f"The expression is {A_k}^{B_k} - {B_k}^{A_k}.")
    print("")
    print("Finding the last digit of each term:")
    print(f"The last digit of {A_k}^{B_k} is {last_digit_A_k_term}.")
    print(f"The last digit of {B_k}^{A_k} is {last_digit_B_k_term}.")
    print("")
    print("The final equation for the last digit is:")
    print(f"({last_digit_A_k_term} - {last_digit_B_k_term}) mod 10 = {final_last_digit}")

solve_last_digit()