def find_last_digit():
    """
    This function calculates the last digit of A_k^(B_k) - B_k^(A_k)
    for a given positive integer k.
    The problem states the last digit is the same for every positive integer k,
    so we can choose any value for k, for instance, k=1.
    """
    k = 1

    # Define A_k and B_k for the chosen k
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    # The expression is A_k^B_k - B_k^A_k.
    # To find the last digit, we perform calculations modulo 10.

    # Calculate the last digit of the first term, A_k^B_k.
    # The pow(base, exponent, modulus) function is highly efficient for this.
    # It calculates (base^exponent) % modulus.
    last_digit_term1 = pow(A_k, B_k, 10)

    # Calculate the last digit of the second term, B_k^A_k.
    last_digit_term2 = pow(B_k, A_k, 10)

    # Calculate the last digit of the final expression.
    # Adding 10 before taking the aodulo ensures the result is positive,
    # which is good practice although Python's % operator handles negative numbers correctly.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    # Output the components of the final equation as requested.
    print(f"For k={k}, we are evaluating the last digit of {A_k}^{B_k} - {B_k}^{A_k}.")
    print(f"The last digit of the first term ({A_k}^{B_k}) is {last_digit_term1}.")
    print(f"The last digit of the second term ({B_k}^{A_k}) is {last_digit_term2}.")
    print(f"The last digit of the difference is ({last_digit_term1} - {last_digit_term2}) mod 10, which equals {final_last_digit}.")

# Run the function to display the result.
find_last_digit()