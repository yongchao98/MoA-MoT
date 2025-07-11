def solve_last_digit():
    """
    This function calculates the last digit of A_k^B_k - B_k^A_k for a given k,
    and explains the logic. The logic is valid for all positive integers k.
    """
    # We can choose any positive integer for k to demonstrate the result.
    # Let's use k = 2.
    k = 2

    # Define A_k and B_k
    # A_k = (10^(k+1) - 1) / 9
    # B_k = 10^k
    # For k=2, A_2 = 111, B_2 = 100
    a_k = (10**(k + 1) - 1) // 9
    b_k = 10**k

    print(f"Let's demonstrate the solution for k = {k}.")
    print(f"A_{k} = (10^({k}+1) - 1) / 9 = {a_k}")
    print(f"B_{k} = 10^{k} = {b_k}")
    print(f"We need to find the last digit of {a_k}^{b_k} - {b_k}^{a_k}.")
    print("-" * 20)

    # To find the last digit, we compute the value modulo 10.
    # The pow(base, exp, mod) function is efficient for this.
    
    # First term: A_k ^ B_k
    # The last digit of A_k (which is always a sequence of 1s) is 1.
    # Any positive integer power of a number ending in 1 will also end in 1.
    last_digit_term1 = pow(a_k, b_k, 10)
    print(f"The last digit of A_{k}^{{B_{k}}} (i.e., {a_k}^{b_k}) is {last_digit_term1}.")

    # Second term: B_k ^ A_k
    # The last digit of B_k (10^k) is 0 for any positive k.
    # Any positive integer power of a number ending in 0 will also end in 0.
    # Note: A_k is always a positive integer.
    last_digit_term2 = pow(b_k, a_k, 10)
    print(f"The last digit of B_{k}^{{A_{k}}} (i.e., {b_k}^{a_k}) is {last_digit_term2}.")
    print("-" * 20)

    # The final equation for the last digit
    # We use (a - b + 10) % 10 to handle potential negative results in modular arithmetic
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10
    
    print("The final equation for the last digits is:")
    print(f"{last_digit_term1} - {last_digit_term2} = {final_last_digit}")
    
    print(f"\nThus, for any positive integer k, the last digit of A_k^B_k - B_k^A_k is always {final_last_digit}.")

solve_last_digit()