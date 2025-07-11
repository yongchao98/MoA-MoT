def find_last_digit():
    """
    This script finds the last digit of the expression A_k^(B_k) - B_k^(A_k),
    where A_k = (10^(k+1)-1)/9 and B_k = 10^k.
    The logic holds for any positive integer k. We will demonstrate with k=1.
    """
    
    # We choose k=1 as a representative example.
    k = 1
    
    # For any k >= 1, A_k = 11...1, so its last digit is 1.
    # A_1 = 11
    last_digit_A_k = 1
    
    # For any k >= 1, B_k = 10^k, so its last digit is 0.
    # B_1 = 10
    last_digit_B_k = 0
    
    # For the first term, A_k^(B_k):
    # Any number ending in 1, raised to a positive integer power, results in a number ending in 1.
    # Since k>=1, B_k is a positive integer.
    last_digit_term1 = 1
    
    # For the second term, B_k^(A_k):
    # Any multiple of 10, raised to a positive integer power, results in a number ending in 0.
    # Since k>=1, A_k is a positive integer.
    last_digit_term2 = 0
    
    # The last digit of the difference is the difference of the last digits (mod 10).
    final_digit = (last_digit_term1 - last_digit_term2) % 10
    
    print("Let's analyze the last digit of the expression A_k^(B_k) - B_k^(A_k).")
    print(f"The last digit of A_k is always {last_digit_A_k}.")
    print(f"The last digit of B_k is always {last_digit_B_k} (for k>=1).")
    print(f"So, the last digit of A_k^(B_k) is {last_digit_term1}.")
    print(f"And the last digit of B_k^(A_k) is {last_digit_term2}.")
    print("\nThe final equation for the last digit of the expression is:")
    print(f"{last_digit_term1} - {last_digit_term2} = {final_digit}")

find_last_digit()