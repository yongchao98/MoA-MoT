import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension of the Dali-esque piano.

    The dimension D is found using a box-counting argument for self-affine fractals,
    which leads to the formula:
    D = log(N * (s_y / s_x)) / log(1 / s_x)
    where:
    - N is the number of self-similar pieces (black keys).
    - s_x is the scaling factor in the x-direction.
    - s_y is the scaling factor in the y-direction.
    """

    # --- Step 1: Define the geometry parameters ---
    keyboard_width = 3.0
    keyboard_height = 1.0

    # There are 7 white keys and 5 black keys in an octave.
    # The 7 white keys are equally sized across the width.
    white_key_width = keyboard_width / 7.0
    
    # Black keys are half as wide as white keys.
    black_key_width = white_key_width / 2.0
    
    # Black keys cover 9/14 of the keyboard height.
    black_key_height = keyboard_height * (9.0 / 14.0)

    # --- Step 2: Identify fractal parameters ---
    # N is the number of black keys in an octave.
    N = 5

    # s_x is the horizontal scaling factor (ratio of black key width to keyboard width)
    s_x = black_key_width / keyboard_width

    # s_y is the vertical scaling factor (ratio of black key height to keyboard height)
    s_y = black_key_height / keyboard_height

    # --- Step 3: Calculate the dimension ---
    # The derived formula from the box-counting argument is:
    # D = log(N * (s_y/s_x)) / log(1/s_x)
    
    # Let's simplify the terms for the final equation display
    anisotropy_ratio = s_y / s_x
    inverse_scaling_factor = 1 / s_x
    numerator_term = N * anisotropy_ratio
    denominator_term = inverse_scaling_factor

    dimension = math.log(numerator_term) / math.log(denominator_term)

    # --- Step 4: Print the results ---
    print("This program calculates the Minkowski-Bouligand dimension of the fractal piano keys.")
    print("\n--- Derivation Steps ---")
    print(f"1. Number of self-similar pieces (N): {N}")
    print(f"2. Horizontal scaling factor (s_x = (3/14)/3): {s_x:.4f}")
    print(f"3. Vertical scaling factor (s_y = 9/14): {s_y:.4f}")
    
    print("\nThe dimension D is given by the formula: D = log(N * s_y/s_x) / log(1/s_x)")
    
    print("\n--- Final Equation ---")
    print(f"N = {N}")
    print(f"s_y/s_x = ({s_y:.4f}) / ({s_x:.4f}) = {anisotropy_ratio:.1f}")
    print(f"1/s_x = 1 / {s_x:.4f} = {denominator_term:.1f}")

    print(f"\nD = log({N} * {anisotropy_ratio:.1f}) / log({denominator_term:.1f})")
    print(f"D = log({numerator_term:.1f}) / log({denominator_term:.1f})")
    
    print(f"\nThe calculated dimension is: {dimension:.4f}")

calculate_fractal_dimension()