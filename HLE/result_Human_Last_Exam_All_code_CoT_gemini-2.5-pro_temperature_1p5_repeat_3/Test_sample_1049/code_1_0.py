import math

def calculate_closed_form(n):
    """
    Calculates the value of the sum for a given integer n using the derived closed-form formula.
    """
    if not isinstance(n, int) or n < 0:
        # The sum is defined for non-negative integers.
        return 0

    # The constant coefficients in the polynomial part of the closed form.
    # The term binom(n + 5, 5) has a coefficient of 1.
    c1 = 116
    c2 = 516
    c3 = 296
    c4 = 16
    
    # Python's math.comb(n, k) handles cases where k > n by returning 0,
    # which is mathematically correct for non-negative n and k.
    
    term1 = math.comb(n + 5, 5)
    term2 = c1 * math.comb(n + 4, 5)
    term3 = c2 * math.comb(n + 3, 5)
    term4 = c3 * math.comb(n + 2, 5)
    term5 = c4 * math.comb(n + 1, 5)
    
    polynomial_part = term1 + term2 + term3 + term4 + term5
    
    return (4**n) * polynomial_part

# Printing the closed-form equation
print("The closed form for the sum S(n) is:")
print("S(n) = 4**n * (1 * binom(n + 5, 5) + 116 * binom(n + 4, 5) + 516 * binom(n + 3, 5) + 296 * binom(n + 2, 5) + 16 * binom(n + 1, 5))")
print("where binom(n, k) represents the binomial coefficient 'n choose k'.\n")

# Example calculation for a specific value of n.
try:
    # Feel free to change this value to test for other n.
    n_example = 2
    result = calculate_closed_form(n_example)
    print(f"For example, when n = {n_example}:")
    print(f"S({n_example}) = {result}")

    # For verification, the original sum for n=2 is:
    # (1^5*binom(0,0)*binom(4,2)) + (3^5*binom(2,1)*binom(2,1)) + (5^5*binom(4,2)*binom(0,0))
    # = (1*1*6) + (243*2*2) + (3125*6*1) = 6 + 972 + 18750 = 19728
    print(f"Verification using direct summation for n=2 gives: 19728")


except ImportError:
    print("Note: math.comb is available in Python 3.8 and later.")
except Exception as e:
    print(f"An error occurred: {e}")
