def solve():
    """
    Calculates the last digit of A_k^(B_k) - B_k^(A_k) for a given k.
    Where A_k = (10^(k+1) - 1) / 9 and B_k = 10^k.
    The problem states the last digit is the same for every positive integer k.
    We will demonstrate this for k=1.
    """
    k = 1  # We choose the smallest positive integer for k.

    # Calculate A_k and B_k
    # For k=1, A_k = (10^2 - 1) / 9 = 11
    Ak = (10**(k + 1) - 1) // 9

    # For k=1, B_k = 10^1 = 10
    Bk = 10**k

    # We need to find the last digit of Ak^Bk - Bk^Ak.
    # This is equivalent to (Ak^Bk - Bk^Ak) mod 10.
    
    # The last digit of Ak^Bk is pow(Ak, Bk, 10).
    # Since the last digit of Ak is 1, the last digit of Ak^Bk is 1.
    last_digit_A_pow_B = pow(Ak, Bk, 10)
    
    # The last digit of Bk^Ak is pow(Bk, Ak, 10).
    # Since the last digit of Bk is 0 (for k>=1), the last digit of Bk^Ak is 0.
    last_digit_B_pow_A = pow(Bk, Ak, 10)
    
    # The last digit of the difference is calculated as:
    # (a - b) mod n = ((a mod n - b mod n) + n) mod n
    final_last_digit = (last_digit_A_pow_B - last_digit_B_pow_A + 10) % 10

    print(f"For k = {k}:")
    print(f"A_k = (10^({k}+1) - 1) / 9 = {Ak}")
    print(f"B_k = 10^{k} = {Bk}")
    print(f"We want to find the last digit of the expression: {Ak}^{Bk} - {Bk}^{Ak}")
    print(f"The last digit of the first term, {Ak}^{Bk}, is {last_digit_A_pow_B}.")
    print(f"The last digit of the second term, {Bk}^{Ak}, is {last_digit_B_pow_A}.")
    print(f"The final equation for the last digit is ({last_digit_A_pow_B} - {last_digit_B_pow_A}) mod 10.")
    print(f"Therefore, the last digit of the entire expression is {final_last_digit}.")

solve()