def solve_last_digit():
    """
    This function determines the last digit of the expression A_k^(B_k) - B_k^(A_k)
    for any positive integer k, where A_k = (10^(k+1)-1)/9 and B_k = 10^k.
    The logic is explained step-by-step.
    """

    print("We want to find the last digit of the expression A_k^B_k - B_k^A_k.")
    print("The last digit of a number N is N mod 10.")
    print("So, we need to compute (A_k^B_k - B_k^A_k) mod 10.")
    print("This is equivalent to ( (A_k^B_k mod 10) - (B_k^A_k mod 10) ) mod 10.\n")

    # Step 1: Analyze the last digit of A_k and B_k
    print("--- Step 1: Find the last digit of A_k and B_k ---")
    print("A_k = (10^(k+1) - 1) / 9 is the number consisting of (k+1) ones. For example, A_1 = 11, A_2 = 111.")
    last_digit_A_k = 1
    print(f"For any positive integer k, the last digit of A_k is always {last_digit_A_k}.")

    print("\nB_k = 10^k. For example, B_1 = 10, B_2 = 100.")
    last_digit_B_k = 0
    print(f"For any positive integer k, B_k is a multiple of 10, so its last digit is always {last_digit_B_k}.\n")

    # Step 2: Find the last digit of the first term, A_k^B_k
    print("--- Step 2: Find the last digit of A_k^B_k ---")
    print("The last digit of A_k^B_k depends on the last digit of A_k.")
    print("Any positive integer power of a number ending in 1 will also end in 1.")
    print(f"Since the last digit of A_k is {last_digit_A_k}, and B_k is a positive integer, the last digit of A_k^B_k must be 1.")
    last_digit_term1 = 1
    print(f"Last digit of A_k^B_k = {last_digit_term1}.\n")

    # Step 3: Find the last digit of the second term, B_k^A_k
    print("--- Step 3: Find the last digit of B_k^A_k ---")
    print("The last digit of B_k^A_k depends on the last digit of B_k.")
    print("Any positive integer power of a number ending in 0 will also end in 0 (since the number is a multiple of 10).")
    print(f"Since the last digit of B_k is {last_digit_B_k}, and A_k is a positive integer, the last digit of B_k^A_k must be 0.")
    last_digit_term2 = 0
    print(f"Last digit of B_k^A_k = {last_digit_term2}.\n")

    # Step 4: Calculate the final result
    print("--- Step 4: Calculate the final result ---")
    print("We need to find the last digit of the difference.")
    print(f"The calculation for the last digit is: (last_digit(A_k^B_k) - last_digit(B_k^A_k)) mod 10")
    
    # We use (a - b + 10) % 10 to handle potential negative results correctly.
    final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10
    
    print(f"The equation for the last digits is: {last_digit_term1} - {last_digit_term2}")
    print(f"The result is {final_last_digit}.")
    print("\nThis last digit is the same for every positive integer k.")

if __name__ == "__main__":
    solve_last_digit()