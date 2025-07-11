import math
from fractions import Fraction

def solve_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the described fractal piano.
    """
    # Step 1: Define the initial parameters using fractions for precision.
    N = 5  # Number of black keys (self-similar copies)
    W = Fraction(3)  # Width of the main keyboard
    H = Fraction(1)  # Height of the main keyboard

    # Step 2: Calculate the dimensions of the black keys.
    # The total width is composed of 7 white keys.
    w_white = W / 7
    # Black keys are half as wide as white keys.
    w_black = Fraction(1, 2) * w_white
    # Black keys have 9/14 of the total height.
    h_black = Fraction(9, 14) * H

    # Step 3: Calculate the anisotropic scaling factors.
    r_x = w_black / W
    r_y = h_black / H

    print("This is a self-affine fractal. We will calculate its box-counting dimension.")
    print("\nStep 1: Determine the parameters of the transformation.")
    print(f"Number of self-similar copies, N = {N}")
    print(f"Scaling factor in x-direction, r_x = ( (1/2) * (3/7) ) / 3 = {r_x}")
    print(f"Scaling factor in y-direction, r_y = (9/14) * 1 = {r_y}")

    # Step 4: Use the formula for the box-counting dimension of this type of self-affine set.
    # Since r_x (1/14) is smaller than r_y (9/14), the formula is:
    print("\nStep 2: State the dimension formula for this self-affine set.")
    print("D = log(N * r_y / r_x) / log(1 / r_x)")

    # Step 5: Substitute the fractional values to simplify the terms.
    term_numerator_frac = N * r_y / r_x
    term_denominator_frac = 1 / r_x

    print("\nStep 3: Simplify the terms inside the logarithms.")
    print(f"Term in numerator's log: N * r_y / r_x = {N} * ({r_y}) / ({r_x}) = {term_numerator_frac}")
    print(f"Term in denominator's log: 1 / r_x = 1 / ({r_x}) = {term_denominator_frac}")

    # Convert to float for the final calculation
    term_num = float(term_numerator_frac)
    term_den = float(term_denominator_frac)

    numerator = math.log(term_num)
    denominator = math.log(term_den)
    dimension = numerator / denominator

    # Step 6: Write the simplified final equation and calculate the result.
    print("\nStep 4: The final equation and result are:")
    print(f"D = log({int(term_num)}) / log({int(term_den)})")
    print(f"D = {numerator:.4f} / {denominator:.4f}")
    print(f"D = {dimension:.4f}")

solve_fractal_dimension()