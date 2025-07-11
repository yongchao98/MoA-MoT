import mpmath
import math

def find_machin_coeffs():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given
    Machin-like formula for pi using the PSLQ algorithm.
    """
    # Set a high precision for the calculations, which is essential for PSLQ.
    mpmath.mp.dps = 150

    # The denominators for the arctan terms in the equation.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of numbers for which to find an integer relation.
    # The equation is of the form: sum(c_i * term_i) - n * (pi/4) = 0
    # So the terms are the arctan values and pi/4.
    terms = [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]
    terms.append(mpmath.pi)

    # Use PSLQ to find a list of integer coefficients C such that C . terms = 0.
    # We add pi to the list of terms because Machin-like relations can sometimes
    # sum to a multiple of pi, not just zero. PSLQ will find the correct
    # linear combination.
    # The relation sought is sum(c_i*atan) - n*pi/4 = k*pi.
    # -> sum(c_i*atan) - (n/4 + k)*pi = 0
    # So we should use `pi` as a basis vector.
    
    # After discovering that simple PSLQ on [atan..., pi/4] may not work if the
    # sum is a multiple of pi, a more robust approach is to include pi
    # itself in the basis vector for PSLQ.
    
    # V = [atan(1/122), ..., atan(1/12943), pi/4, pi]
    extended_terms = [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]
    extended_terms.append(mpmath.pi / 4)
    extended_terms.append(mpmath.pi)
    
    coeffs_ext = mpmath.pslq(extended_terms)
    
    # The found relation is:
    # C[0]t_0 + ... + C[5]t_5 + C[6]*(pi/4) + C[7]*pi = 0
    # We want to match this to: c_1*t_1 + ... + c_6*t_6 = n*pi/4
    # c_1*t_1 + ... - n*pi/4 = 0
    # Our found relation gives: C[0]t_0 + ... = -C[6]*pi/4 - C[7]*pi = (-C[6] - 4*C[7])*pi/4
    # Thus, n = -C[6] - 4*C[7] and c_i = C[i]

    int_coeffs = [int(round(c)) for c in coeffs_ext]

    c = int_coeffs[0:6]
    n_val = -int_coeffs[6] - 4 * int_coeffs[7]

    # Normalize for smallest positive n.
    if n_val < 0:
        n_val = -n_val
        c = [-x for x in c]

    # Make the relation primitive by dividing by GCD.
    def gcd_list(numbers):
        if not numbers:
            return 1
        result = abs(numbers[0])
        for i in range(1, len(numbers)):
            result = math.gcd(result, abs(numbers[i]))
        return result

    all_coeffs_for_gcd = c + [n_val]
    common_divisor = gcd_list(all_coeffs_for_gcd)
    
    if common_divisor > 1:
        n_val //= common_divisor
        c = [x // common_divisor for x in c]

    # Print the equation with the found coefficients.
    print(f"The discovered relation is:")
    
    # Display n * pi / 4
    lhs = f"{n_val}*pi/4"
    
    # Build the RHS string
    rhs_parts = []
    first = True
    for i in range(len(c)):
        if c[i] == 0:
            continue
        
        # Add sign
        if not first:
            if c[i] > 0:
                rhs_parts.append("+ ")
            else:
                rhs_parts.append("- ")
        elif c[i] < 0:
            rhs_parts.append("-")
            
        # Add coefficient (if not 1)
        coeff_abs = abs(c[i])
        if coeff_abs != 1:
            rhs_parts.append(f"{coeff_abs}*")
            
        rhs_parts.append(f"arctan(1/{denominators[i]}) ")
        first = False
        
    print(f"{lhs} = {''.join(rhs_parts)}")

    # Final Answer in requested format
    final_coeffs = [n_val] + c
    print("\nCoefficients (n, c1, ..., c6):")
    print(",".join(map(str, final_coeffs)))
    
    # Return final answer per instruction format.
    # Format: n,c1,c2,c3,c4,c5,c6
    final_answer_str = ",".join(map(str, final_coeffs))
    print(f"\n<<< {final_answer_str} >>>")


find_machin_coeffs()
