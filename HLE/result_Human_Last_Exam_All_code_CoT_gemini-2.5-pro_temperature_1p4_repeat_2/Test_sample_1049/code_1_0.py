import math

def print_closed_form():
    """
    This function prints the derived closed-form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    where C(n, k) is the binomial coefficient.
    
    The script will output the formula and its constituent numerical parts.
    """

    # Coefficients of the quartic polynomial Q(n) in the closed-form expression.
    c4 = 63
    c3 = 182
    c2 = 173
    c1 = 62
    c0 = 8

    # The closed form is: S_n = 2**(2n - 3) * (n + 1) * Q(n)
    
    # Constructing the string for the final equation to present to the user.
    quartic_poly_str = f"{c4}*n**4 + {c3}*n**3 + {c2}*n**2 + {c1}*n + {c0}"
    final_eq_str = f"S_n = 2**(2*n - 3) * (n + 1) * ({quartic_poly_str})"

    print("The closed form for the sum is given by the following equation:")
    print(final_eq_str)

    print("\nTo satisfy the request to output each number in the final equation, here is a breakdown:")
    
    print("\nThe formula can be written as S_n = (Factor 1) * (Factor 2) * (Factor 3)")
    print("Factor 1 (a power of 2): 2**(2*n - 3)")
    print("Factor 2 (a linear term): (n + 1)")
    print(f"Factor 3 (a quartic polynomial in n): {quartic_poly_str}")

    print("\nThe coefficients of the quartic polynomial are:")
    print(f"  - Coefficient of n**4: {c4}")
    print(f"  - Coefficient of n**3: {c3}")
    print(f"  - Coefficient of n**2: {c2}")
    print(f"  - Coefficient of n**1: {c1}")
    print(f"  - The constant term: {c0}")

if __name__ == '__main__':
    print_closed_form()