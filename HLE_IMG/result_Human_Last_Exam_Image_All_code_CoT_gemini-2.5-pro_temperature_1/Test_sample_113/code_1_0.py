import math
from fractions import Fraction

def calculate_fractal_circle_area():
    """
    This script calculates the limit of the total area of circles in a fractal pattern.
    The problem is solved by identifying and summing an infinite geometric series.
    """

    # Step 1: Analyze the initial figure (R1) and find the area of the first circle.
    print("--- Step 1: Area of the first circle ---")
    w0 = 6
    h0 = 8
    print(f"The initial rectangle has a width w0 = {w0} and a height h0 = {h0}.")

    # The radius of the first circle (r1) is one-third of the initial width.
    r1 = Fraction(w0, 3)
    print(f"The radius of the first circle is r1 = w0 / 3 = {w0}/3 = {r1}.")

    # The area of the first circle (A1) is pi * r1^2. We'll track the coefficient of pi.
    A1_coeff = r1**2
    print(f"The area of this circle is A1 = pi * r1^2 = pi * ({r1})^2 = {A1_coeff}*pi.")
    print("-" * 40)

    # Step 2: Determine the geometric scaling factor between generations.
    print("--- Step 2: Scaling factor between generations ---")
    diagonal_length = math.sqrt(w0**2 + h0**2)
    print(f"The length of the diagonal of the initial rectangle is sqrt({w0}^2 + {h0}^2) = {int(diagonal_length)}.")

    # The linear scaling factor (k) for the dimensions of the rectangles is derived from the geometry.
    # A new rectangle's width w_new is related to the previous width w_old by k = w_new/w_old.
    # The derivation gives k = 1/2 - w0 / (3 * diagonal_length).
    k_lin = Fraction(1, 2) - Fraction(w0, 3 * int(diagonal_length))
    print(f"The linear scaling factor for dimensions is k = 1/2 - {w0}/(3*{int(diagonal_length)}) = {k_lin}.")
    print("-" * 40)

    # Step 3: Formulate the total area as a geometric series.
    print("--- Step 3: Setting up the geometric series ---")
    # At each step, 4 new circles are created from each previous one.
    # The area of each new circle is scaled by k^2.
    # So, the common ratio of the total area series is r = 4 * k^2.
    ratio = 4 * k_lin**2
    print(f"The common ratio of the geometric series for total area is r = 4 * k^2 = 4 * ({k_lin})^2 = {ratio}.")
    print("-" * 40)

    # Step 4: Sum the infinite geometric series to find the total area.
    print("--- Step 4: Summing the series ---")
    # The total area S is the sum of the infinite series S = a / (1 - r),
    # where the first term is a = A1 and the common ratio is r.
    print(f"The total area S = A1 / (1 - r).")
    print(f"The first term is a = {A1_coeff}*pi.")
    print(f"The common ratio is r = {ratio}.")

    # Perform the final calculation. S_coeff = A1_coeff / (1 - ratio).
    S_coeff = A1_coeff / (1 - ratio)

    # Print the final equation with all numbers.
    print("\nFinal Equation Calculation:")
    # S = (A1) / (1 - r)
    print(f"S = ({A1_coeff.numerator}/{A1_coeff.denominator} * pi) / (1 - {ratio.numerator}/{ratio.denominator})")
    
    # S = (A1) / ( (r_den - r_num)/r_den )
    denominator_val = 1 - ratio
    print(f"S = ({A1_coeff.numerator}/{A1_coeff.denominator} * pi) / ({denominator_val.numerator}/{denominator_val.denominator})")
    
    # S = (A1_num * den_den * pi) / (A1_den * den_num)
    unsimplified_num = A1_coeff.numerator * denominator_val.denominator
    unsimplified_den = A1_coeff.denominator * denominator_val.numerator
    print(f"S = ({unsimplified_num}*pi) / {unsimplified_den}")

    # Simplify the final fraction.
    simplified_num = S_coeff.numerator
    simplified_den = S_coeff.denominator
    print(f"\nAfter simplification, the final result is:")
    print(f"S = ({simplified_num}*pi) / {simplified_den}")


calculate_fractal_circle_area()