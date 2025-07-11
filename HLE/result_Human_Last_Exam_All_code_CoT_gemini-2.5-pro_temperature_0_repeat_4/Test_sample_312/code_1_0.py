import fractions

def solve_frostman_decay():
    """
    Calculates the smallest possible decay exponent 'c' based on the provided problem.
    """
    # Parameters from the problem statement
    n = 2
    alpha_num = 8
    alpha_den = 5
    alpha = fractions.Fraction(alpha_num, alpha_den)

    # The critical exponent c is determined by a well-known formula in Fourier analysis
    # related to the spherical averaging of Fourier transforms of measures.
    # The formula is: c = min(-(n - alpha), -(n - 1) / 2)

    # Calculate the two components of the formula using fractions for precision
    term1 = -(n - alpha)
    term2 = fractions.Fraction(-(n - 1), 2)

    # Determine the minimum to find c
    c = min(term1, term2)

    # --- Output the step-by-step calculation ---

    print("The smallest possible value for c is determined by the formula from Fourier analysis:")
    print("c = min(-(n - alpha), -(n - 1) / 2)\n")

    print(f"Given parameters are n = {n} and alpha = {alpha.numerator}/{alpha.denominator}.\n")

    print("Plugging the values into the formula:")
    # Print the equation with the numbers substituted
    print(f"c = min(-({n} - {alpha.numerator}/{alpha.denominator}), -({n} - 1)/2)")
    print(f"c = min(-({n * alpha.denominator - alpha.numerator}/{alpha.denominator}), {term2.numerator}/{term2.denominator})")
    print(f"c = min({term1.numerator}/{term1.denominator}, {term2.numerator}/{term2.denominator})")
    
    # Convert fractions to float for comparison in the print statement
    term1_float = float(term1)
    term2_float = float(term2)
    c_float = float(c)
    
    print(f"c = min({term1_float}, {term2_float})\n")

    print(f"Comparing the two values, we find the minimum is {c_float}.")
    print(f"Therefore, the smallest possible value for c is {c.numerator}/{c.denominator}.")

solve_frostman_decay()