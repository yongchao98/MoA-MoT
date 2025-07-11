import mpmath

def find_machin_like_coefficients():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given
    Machin-like formula for pi.
    
    The equation is:
    n*pi/4 = c1*arctan(1/122) + c2*arctan(1/239) + c3*arctan(1/682) + 
             c4*arctan(1/1252) + c5*arctan(1/2855) + c6*arctan(1/12943)

    This is solved as an integer relation problem using the PSLQ algorithm.
    """
    # Set a high precision for the calculations, as required by PSLQ.
    mpmath.mp.dps = 100

    # These are the denominators in the arctan terms.
    x_vals = [122, 239, 682, 1252, 2855, 12943]

    # We are looking for an integer relation for the equation:
    # n*(pi/4) - c1*arctan(1/122) - c2*arctan(1/239) - ... = 0
    # We construct a vector of the constant values.
    print("Constructing the vector of values for PSLQ...")
    v = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / x) for x in x_vals]

    # Run the PSLQ algorithm to find a vector of integer coefficients [k0, k1, ...]
    # such that k0*v[0] + k1*v[1] + ... = 0.
    print("Running PSLQ algorithm... (this may take a moment)")
    coeffs = mpmath.pslq(v)
    print("PSLQ has found an integer relation.")

    # From the relation found by PSLQ, we extract our coefficients:
    # n = k0
    # c_i = -k_i for i > 0
    n = coeffs[0]
    c = [-val for val in coeffs[1:]]

    # The problem requires the smallest *positive* n. If PSLQ returns a
    # negative n, we can multiply the entire relation by -1 to make n positive.
    if n < 0:
        n = -n
        c = [-val for val in c]

    # Display the final discovered equation.
    print("\n--------------------------------------------------")
    print("The discovered relation is:")
    
    terms = []
    for i in range(len(c)):
        if c[i] != 0:
            # Format coefficient, hiding it if it's 1 or -1.
            if abs(c[i]) == 1:
                coeff_str = "" if c[i] > 0 else "-"
            else:
                coeff_str = f"{abs(c[i])}*"

            sign = " + " if c[i] > 0 else " - "
            
            terms.append((sign, f"{coeff_str}arctan(1/{x_vals[i]})"))

    # Nicely format the equation string
    equation_str = f"{n} * (pi/4) ="
    # Handle the first term separately to avoid a leading "+"
    equation_str += f" {terms[0][1]}" if terms[0][0] == " + " else f" -{terms[0][1]}"
    for sign, term in terms[1:]:
        equation_str += f"{sign}{term}"
        
    print(equation_str)
    
    # Print the values of n and c_i as the final answer.
    print("\n--------------------------------------------------")
    print("The integer coefficients are:")
    print(f"n = {n}")
    for i, val in enumerate(c):
        print(f"c_{i+1} = {val}")

    # Final answer in the specified format
    final_answer = f"{n},{','.join(map(str, c))}"
    print("\nFinal Answer:")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    find_machin_like_coefficients()