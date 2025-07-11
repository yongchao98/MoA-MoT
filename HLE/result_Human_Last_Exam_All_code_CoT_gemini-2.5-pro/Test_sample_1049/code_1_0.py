import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_closed_form(n):
    """
    Calculates the sum using the derived closed-form formula.
    S_n = 4^n * (C(n+5,5) + 116*C(n+4,5) + 516*C(n+3,5) + 296*C(n+2,5) + 16*C(n+1,5))
    """
    if n < 0:
        return 0
        
    term1 = combinations(n + 5, 5)
    term2 = 116 * combinations(n + 4, 5)
    term3 = 516 * combinations(n + 3, 5)
    term4 = 296 * combinations(n + 2, 5)
    term5 = 16 * combinations(n + 1, 5)
    
    polynomial_part = term1 + term2 + term3 + term4 + term5
    
    return (4**n) * polynomial_part

def print_formula():
    """
    Prints the derived closed-form formula for the sum.
    """
    print("The closed form for the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k,n-k) is:")
    print("\nS_n = 4^n * P(n)\n")
    print("where P(n) is a polynomial in n, which can be expressed as:")
    print("P(n) = C(n+5, 5) + 116*C(n+4, 5) + 516*C(n+3, 5) + 296*C(n+2, 5) + 16*C(n+1, 5)")
    print("\nHere, C(n, k) denotes the binomial coefficient 'n choose k'.")


if __name__ == '__main__':
    # Print the formula
    print_formula()
    
    # Example calculation for a specific n
    n_value = 5
    result = calculate_closed_form(n_value)
    
    print(f"\nFor n = {n_value}, the value of the sum is:")
    print(result)

    # Verification for n=1
    # S_1 = (1^5 * C(0,0)*C(2,1)) + (3^5 * C(2,1)*C(0,0))
    #     = (1 * 1 * 2) + (243 * 2 * 1) = 2 + 486 = 488
    # Using the formula:
    # S_1 = 4^1 * (C(6,5) + 116*C(5,5)) = 4 * (6 + 116) = 4 * 122 = 488
    print("\n--- Verification for n=1 ---")
    print(f"Direct calculation for n=1 gives 488.")
    print(f"Formula calculation for n=1 gives: {calculate_closed_form(1)}")

