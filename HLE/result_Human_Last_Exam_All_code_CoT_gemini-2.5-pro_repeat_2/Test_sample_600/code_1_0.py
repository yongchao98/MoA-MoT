def solve():
    """
    Calculates the last digit of A_k^B_k - B_k^A_k.
    The problem states the last digit is the same for every positive integer k.
    We will demonstrate this for a sample value, k=3.
    """
    k = 3

    # Define A_k and B_k
    # A_k = (10**(k+1) - 1) / 9 is a repunit with k+1 ones.
    A_k = (10**(k + 1) - 1) // 9
    # B_k = 10**k
    B_k = 10**k

    print(f"Let's demonstrate with k = {k}:")
    print(f"A_{k} = (10^({k}+1) - 1) / 9 = {A_k}")
    print(f"B_{k} = 10^{k} = {B_k}")
    print(f"We need the last digit of A_{k}^(B_{k}) - B_{k}^(A_{k}).")
    print("-" * 30)

    # To find the last digit, we compute the value modulo 10.
    # We use python's pow(base, exp, mod) for efficient modular exponentiation.

    # Last digit of the first term: A_k^(B_k)
    last_digit_term1 = pow(A_k, B_k, 10)
    print(f"The last digit of A_{k}^(B_{k}) is ({A_k} ^ {B_k}) mod 10 = {last_digit_term1}")

    # Last digit of the second term: B_k^(A_k)
    # Since B_k is a multiple of 10 and A_k is a positive integer, this will be 0.
    last_digit_term2 = pow(B_k, A_k, 10)
    print(f"The last digit of B_{k}^(A_{k}) is ({B_k} ^ {A_k}) mod 10 = {last_digit_term2}")
    print("-" * 30)

    # The last digit of the difference is (last_digit_term1 - last_digit_term2) mod 10.
    final_last_digit = (last_digit_term1 - last_digit_term2) % 10

    print("The final equation for the last digit is:")
    print(f"({last_digit_term1} - {last_digit_term2}) mod 10 = {final_last_digit}")
    print(f"\nThus, the constant last digit is {final_last_digit}.")

solve()