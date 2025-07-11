# This script calculates the last digit of A_k^(B_k) - B_k^(A_k) for k=1, 2, and 3
# to verify the result that the last digit is constant.

for k in range(1, 4):
    # Define A_k and B_k for the current value of k.
    # A_k is a sequence of k+1 ones, e.g., 11, 111, etc.
    A_k = (10**(k + 1) - 1) // 9
    # B_k is 1 followed by k zeros, e.g., 10, 100, etc.
    B_k = 10**k

    # Python's pow(base, exp, mod) is very efficient for finding the last digit of large powers.
    # It calculates (base^exp) % mod.
    
    # Calculate the last digit of the first term, A_k^(B_k).
    last_digit_term1 = pow(A_k, B_k, 10)

    # Calculate the last digit of the second term, B_k^(A_k).
    # Since k>=1, B_k is a multiple of 10 and A_k is a positive exponent.
    last_digit_term2 = pow(B_k, A_k, 10)

    # The last digit of (Term1 - Term2) is (last_digit_term1 - last_digit_term2) mod 10.
    # We add 10 before taking the modulo to handle potential negative results correctly.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    # Print the equation with its components as requested, and the final result.
    print(f"For k={k}, we find the last digit of the expression {A_k}^{B_k} - {B_k}^{A_k}")
    print(f"The result is {final_last_digit}.")
    # Add a separator for readability, except for the last line.
    if k < 3:
        print("---")