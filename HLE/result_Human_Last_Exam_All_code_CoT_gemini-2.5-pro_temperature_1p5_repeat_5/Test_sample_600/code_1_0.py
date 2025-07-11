def solve_last_digit():
    """
    Calculates the last digit of A_k^B_k - B_k^A_k for a given k.
    """
    # We can choose any positive integer for k. Let's use k=2 as an example.
    # The logic is the same for all positive integers k.
    k = 2

    # A_k is a sequence of k+1 ones.
    Ak = int('1' * (k + 1))
    # B_k is 10 to the power of k.
    Bk = 10**k

    print(f"The problem states that for every positive integer k, the last digit of A_k^B_k - B_k^A_k is the same.")
    print(f"Let's demonstrate this with k = {k}.")
    print(f"\nFirst, we define A_k and B_k:")
    print(f"A_k = (10^(k+1) - 1) / 9")
    print(f"B_k = 10^k")
    print(f"\nFor k = {k}:")
    print(f"A_{k} = {Ak}")
    print(f"B_{k} = {Bk}")
    print("\nThe expression is A_k^B_k - B_k^A_k, which for our k becomes:")
    print(f"{Ak}^{Bk} - {Bk}^{Ak}")
    print("\nTo find the last digit, we can use modular arithmetic (find the remainder when divided by 10).")

    # Step 1: Find the last digit of the first term, Ak^Bk
    last_digit_Ak = Ak % 10
    # Any positive integer power of a number ending in 1 also ends in 1.
    last_digit_term1 = 1

    print(f"\nStep 1: Analyze the last digit of the first term, {Ak}^{Bk}")
    print(f"The last digit of A_{k} ({Ak}) is {last_digit_Ak}.")
    print(f"Any positive integer power of a number ending in 1 will also end in 1.")
    print(f"Since B_{k} ({Bk}) is a positive integer, the last digit of {Ak}^{Bk} is {last_digit_term1}.")

    # Step 2: Find the last digit of the second term, Bk^Ak
    last_digit_Bk = Bk % 10
    # For k>=1, Bk is a multiple of 10. For any positive integer power Ak, Bk^Ak will end in 0.
    last_digit_term2 = 0

    print(f"\nStep 2: Analyze the last digit of the second term, {Bk}^{Ak}")
    print(f"The last digit of B_{k} ({Bk}) is {last_digit_Bk}.")
    print(f"Any positive integer power of a number ending in 0 will also end in 0.")
    print(f"Since A_{k} ({Ak}) is a positive integer, the last digit of {Bk}^{Ak} is {last_digit_term2}.")
    
    # Step 3: Calculate the last digit of the difference
    final_last_digit = (last_digit_term1 - last_digit_term2) % 10

    print(f"\nStep 3: Calculate the last digit of the difference")
    print(f"The last digit of ({Ak}^{Bk} - {Bk}^{Ak}) is the last digit of (a number ending in {last_digit_term1}) - (a number ending in {last_digit_term2}).")
    print(f"This is calculated as ({last_digit_term1} - {last_digit_term2}) mod 10.")

    print(f"\nThe result is {final_last_digit}.")
    print(f"\nThis confirms that for any positive integer k, the last digit of the expression is always {final_last_digit}.")

solve_last_digit()