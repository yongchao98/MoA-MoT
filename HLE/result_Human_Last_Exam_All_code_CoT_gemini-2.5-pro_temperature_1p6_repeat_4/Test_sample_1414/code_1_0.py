import mpmath
import math

def solve_machin_like_formula():
    """
    This function solves for the integer coefficients n, c_1, ..., c_6
    in the given Machin-like formula for pi.
    
    The equation is:
    n*pi/4 = c1*arctan(1/122) + c2*arctan(1/239) + c3*arctan(1/682) +
             c4*arctan(1/1252) + c5*arctan(1/2855) + c6*arctan(1/12943)

    The method used is the PSLQ algorithm, which is designed to find integer
    relations between real numbers. We are looking for an integer relation for
    the given arctan values and pi/4.
    """
    
    # Set high precision for PSLQ to work reliably.
    mpmath.mp.dps = 200

    x_values = [122, 239, 682, 1252, 2855, 12943]
    
    # The list of values for which we seek an integer relation.
    terms = [mpmath.atan(mpmath.mpf(1)/x) for x in x_values]
    terms.append(mpmath.pi / 4)

    # Use mpmath's PSLQ implementation (identify) to find integer coefficients.
    # The relation found will be c1*t1 + c2*t2 + ... + c7*t7 = 0,
    # where t7 is pi/4. We are solving for c1*t1 + ... = n*pi/4, so n = -c7.
    coeffs = mpmath.identify(terms, tol=1e-180)
    
    if not coeffs:
        print("Could not find an integer relation. Try increasing the precision.")
        return

    n = -coeffs[-1]
    cs = list(coeffs[:-1])
    
    # The problem asks for the smallest positive n, so if n is negative,
    # we flip the signs of all coefficients.
    if n < 0:
        n = -n
        cs = [-c for c in cs]

    # To get the solution with the smallest n, we divide all coefficients
    # by their greatest common divisor.
    all_coeffs = cs + [n]
    
    common_divisor = n
    for c in cs:
        common_divisor = math.gcd(common_divisor, c)
        
    n = n // common_divisor
    c1, c2, c3, c4, c5, c6 = [c // common_divisor for c in cs]
    
    print("The discovered identity is:")
    print(f"{n} * \u03C0/4 =", end="")
    
    final_coeffs = [c1, c2, c3, c4, c5, c6]
    is_first_term = True
    for i in range(len(x_values)):
        if final_coeffs[i] != 0:
            val = abs(final_coeffs[i])
            # Set sign, but handle first term differently
            if final_coeffs[i] > 0:
                sign = " + "
                if is_first_term:
                    sign = " "
            else:
                sign = " - "

            print(f"{sign}{val}*arctan(1/{x_values[i]})", end="")
            is_first_term = False
    print("\n")

    print("The constants are:")
    print(f"n,c1,c2,c3,c4,c5,c6 = {n},{c1},{c2},{c3},{c4},{c5},{c6}")

if __name__ == '__main__':
    solve_machin_like_formula()