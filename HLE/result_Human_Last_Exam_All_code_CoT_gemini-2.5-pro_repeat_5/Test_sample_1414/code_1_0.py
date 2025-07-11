import mpmath

def solve_pi_formula():
    """
    Finds the integer coefficients for a Machin-like formula for pi using the PSLQ algorithm.
    The formula is of the form:
    n * pi / 4 = c1*arctan(1/122) + c2*arctan(1/239) + c3*arctan(1/682) + 
                   c4*arctan(1/1252) + c5*arctan(1/2855) + c6*arctan(1/12943)
    """

    # Set the precision for the calculations. PSLQ requires high precision.
    mpmath.mp.dps = 100

    # The denominators for the arctan arguments
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of high-precision real numbers for the PSLQ algorithm.
    # The vector is [arctan(1/x1), arctan(1/x2), ..., arctan(1/x6), pi/4]
    x = [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]
    x.append(mpmath.pi / 4)

    # Find the integer relation using PSLQ.
    # This finds a list of integers 'coeffs' such that:
    # coeffs[0]*x[0] + ... + coeffs[6]*x[6] = 0
    coeffs = mpmath.pslq(x)

    # The problem is stated as:
    # c1*x[0] + ... + c6*x[5] = n * (pi/4)
    # Rearranging this gives:
    # c1*x[0] + ... + c6*x[5] - n * (pi/4) = 0
    # Our relation is:
    # coeffs[0]*x[0] + ... + coeffs[5]*x[5] + coeffs[6]*x[6] = 0
    # where x[6] is pi/4.
    # So, we can identify c_k = coeffs[k] for k=0..5 and -n = coeffs[6].
    
    c = coeffs[:-1]
    n_coeff = coeffs[-1]

    # Therefore, n = -n_coeff
    n = -n_coeff

    # The problem asks for the smallest positive n.
    # If PSLQ returns a relation where n would be negative, we can multiply
    # the entire relation by -1 to make n positive. This just flips the signs
    # of all coefficients.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    # Unpack coefficients for printing
    c1, c2, c3, c4, c5, c6 = c

    # Print the full equation as requested, showing each number
    print("The solved equation is:")
    print(f"{n}*pi/4 = {c1}*arctan(1/{denominators[0]}) + {c2}*arctan(1/{denominators[1]}) + {c3}*arctan(1/{denominators[2]}) + {c4}*arctan(1/{denominators[3]}) + {c5}*arctan(1/{denominators[4]}) + {c6}*arctan(1/{denominators[5]})")

    print("\nThe constants are:")
    print(f"n = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")
    
    # Prepare the final answer string for the <<<...>>> format
    final_answer_list = [n, c1, c2, c3, c4, c5, c6]
    return ",".join(map(str, final_answer_list))

# Execute the function to find the solution and get the final answer string.
final_answer = solve_pi_formula()

# The final answer must be wrapped in <<<>>>
print(f"\n<<<{final_answer}>>>")
