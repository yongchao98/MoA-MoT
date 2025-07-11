import mpmath
import math

def gcd_list(numbers):
    """
    Computes the greatest common divisor of a list of integers.
    """
    if not numbers:
        return 0
    # The list might contain zeros, filter them out as gcd(k,0)=|k|.
    # math.gcd requires non-negative integers as of Python 3.9+.
    non_zero_numbers = [abs(n) for n in numbers if n != 0]
    if not non_zero_numbers:
        return 1 # GCD of empty set/all zeros is ill-defined, 1 is a safe identity.

    result = non_zero_numbers[0]
    for i in range(1, len(non_zero_numbers)):
        result = math.gcd(result, non_zero_numbers[i])
    return result

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula.
    """
    # Set the precision for mpmath. PSLQ needs high precision.
    mpmath.mp.dps = 100

    # The denominators from the problem statement
    y = [122, 239, 682, 1252, 2855, 12943]

    # The equation is: sum(c_i * atan(1/y_i)) = n * pi/4
    # We reformulate it as: sum(c_i * atan(1/y_i)) - n * pi/4 = 0
    # The vector for PSLQ is [atan(1/y_1), ..., atan(1/y_6), -pi/4]
    # PSLQ will find coefficients [c_1, ..., c_6, n]
    
    vals = [mpmath.atan(mpmath.mpf(1)/v) for v in y]
    vals.append(-mpmath.pi / 4)

    # Run the PSLQ algorithm to find the integer relation
    raw_coeffs = mpmath.pslq(vals)

    # Extract the coefficients c_i and n
    c_coeffs = raw_coeffs[:-1]
    n_coeff = raw_coeffs[-1]

    # The problem asks for the smallest positive n.
    # If n is negative, flip the signs of all coefficients.
    if n_coeff < 0:
        n_coeff = -n_coeff
        c_coeffs = [-c for c in c_coeffs]

    # To get the smallest n, we find the primitive solution
    # by dividing all coefficients by their greatest common divisor.
    all_coeffs = list(c_coeffs) + [n_coeff]
    common_divisor = gcd_list(all_coeffs)
    
    n_final = n_coeff // common_divisor
    c_final = [c // common_divisor for c in c_coeffs]

    # Print the full equation for verification
    print("The discovered identity is:")
    equation_parts = []
    for i, c in enumerate(c_final):
        if c != 0:
            sign = "+" if c > 0 else "-"
            # handle the first term not having a sign
            if not equation_parts:
                sign = "" if c > 0 else "-"
            else:
                 equation_parts.append(f" {sign} ")
            
            abs_c = abs(c)
            term = f"{abs_c}*arctan(1/{y[i]})" if abs_c != 1 else f"arctan(1/{y[i]})"
            equation_parts.append(term)
            
    print(f"{n_final}*pi/4 = {''.join(equation_parts)}")
    
    # Print the coefficients as requested
    print("\nThe solution is:")
    print(f"n = {n_final}")
    for i, c in enumerate(c_final):
        print(f"c_{i+1} = {c}")

    # Provide the final answer in the requested format for parsing.
    print("\nFinal Answer String:")
    final_answer_string = f"{n_final}," + ",".join(map(str, c_final))
    print(f"n,c_1,c_2,c_3,c_4,c_5,c_6 = {final_answer_string}")

if __name__ == '__main__':
    solve_arctan_relation()
