def solve():
    """
    This function calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given k.
    The logic holds for any positive integer k, so we use k=2 as an example.
    """
    k = 2

    # Define A_k and B_k
    # A_k = (10^(k+1) - 1) / 9
    # B_k = 10^k
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k

    # We need the last digit of A_k^(B_k) - B_k^(A_k).
    # This is equivalent to (A_k^(B_k) - B_k^(A_k)) mod 10.
    
    # Python's pow(base, exp, mod) is efficient for modular exponentiation.
    
    # Last digit of the first term: A_k^(B_k)
    last_digit_term1 = pow(A_k, B_k, 10)
    
    # Last digit of the second term: B_k^(A_k)
    last_digit_term2 = pow(B_k, A_k, 10)
    
    # The final result is the last digit of the difference.
    # We add 10 before taking the modulo to handle potential negative results (e.g., (1-2)%10).
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

    print(f"Let's verify the result for k = {k}:")
    print(f"A_{k} = {A_k}")
    print(f"B_{k} = {B_k}")
    print(f"The last digit of A_{k}^B_{k} (i.e., {A_k}^{B_k}) is {last_digit_term1}.")
    print(f"The last digit of B_{k}^A_{k} (i.e., {B_k}^{A_k}) is {last_digit_term2}.")
    
    # The final equation showing the numbers involved in the last step
    print("\nThe final equation for the last digit is:")
    print(f"({last_digit_term1} - {last_digit_term2}) mod 10 = {final_last_digit}")
    
    print("\nSince this is true for any positive integer k, the constant last digit is found.")
    print(f"The final answer is: {final_last_digit}")

solve()