from fractions import Fraction

def find_decay_exponent():
    """
    Calculates the smallest possible value for the decay exponent 'c' for the
    L2 norm of a Frostman measure's Fourier transform on a circle.
    """
    # The problem is to find the sharp exponent 'c' for the decay of the
    # spherically averaged Fourier transform of a measure. The general formula
    # for this exponent is c = -(n - s) / 2, where n is the dimension of the
    # ambient space and s is the dimension of the measure.

    # Parameters from the problem statement
    n = 2  # The space is R^2
    s = Fraction(8, 5)  # The measure is a 8/5-Frostman measure

    print("The formula for the smallest possible exponent c is: c = -(n - s) / 2")
    print(f"Here, the spatial dimension is n = {n}.")
    print(f"The dimension of the Frostman measure is s = {s.numerator}/{s.denominator}.")
    print("-" * 50)

    # Step 1: Calculate the difference (n - s)
    n_minus_s = n - s
    n_frac = Fraction(n)

    # Step 2: Calculate the final exponent c
    c = -n_minus_s / 2

    # Step 3: Print the full calculation, showing each number
    print("The step-by-step calculation is as follows:")
    print(f"c = -({n} - {s.numerator}/{s.denominator}) / 2")
    print(f"  = -({n_frac.numerator}/{n_frac.denominator} - {s.numerator}/{s.denominator}) / 2")
    print(f"  = -({n_minus_s.numerator}/{n_minus_s.denominator}) / 2")
    print(f"  = {c.numerator}/{c.denominator}")
    print("-" * 50)
    
    print(f"The smallest possible value for c is {c.numerator}/{c.denominator}.")

# Run the calculation
find_decay_exponent()