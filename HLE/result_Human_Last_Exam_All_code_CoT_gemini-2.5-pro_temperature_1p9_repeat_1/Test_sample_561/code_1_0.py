import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension for the described fractal piano keys.
    """
    # N: Number of smaller pianos (one for each black key).
    # There are 5 black keys in a standard octave.
    N = 5

    # s_x: The scaling factor in the x-direction (width).
    # The full keyboard is 3 units wide.
    # There are 7 white keys, so each white key is 3/7 units wide.
    # Black keys are half as wide: (3/7)/2 = 3/14 units.
    # The scaling factor s_x is the new width divided by the old width.
    # s_x = (3/14) / 3 = 1/14.
    s_x_inv = 14

    # s_y: The scaling factor in the y-direction (height).
    # The full keyboard is 1 unit high.
    # Black keys have a height of 9/14 units.
    # The scaling factor s_y is the new height divided by the old height.
    # s_y = (9/14) / 1 = 9/14.
    s_y_inv = 14 / 9

    # N_x: The number of distinct horizontal positions (columns).
    # The 5 black keys are at 5 different horizontal positions.
    N_x = 5
    # N_y: The number of distinct vertical positions (rows).
    # The 5 black keys are aligned in the same vertical band.
    N_y = 1
    
    # The dimension 'D' for a self-affine set of this type is given by the formula:
    # D = log(N_y) / log(1/s_y) + log(N_x) / log(1/s_x)
    # Plugging in the values:
    # D = log(1) / log(14/9) + log(5) / log(14)
    # Since log(1) = 0, the first term vanishes.
    # D = 0 + log(5) / log(14)
    
    # We will now calculate the final term.
    numerator = N_x # which is N / N_y = 5 / 1 = 5
    denominator_base = s_x_inv

    dimension = math.log(numerator) / math.log(denominator_base)
    
    print("The Minkowski–Bouligand dimension (D) for this self-affine fractal is calculated using the formula:")
    print("D = log(N_y) / log(1/s_y) + log(N_x) / log(1/s_x)")
    print("\nBased on the problem description:")
    print(f"- Number of recursive pieces (N) = {N}")
    print(f"- Horizontal scaling factor (s_x) = 1/{s_x_inv}")
    print(f"- Vertical scaling factor (s_y) = 9/{s_y_inv * 9}")
    print(f"- Number of distinct horizontal placements (N_x) = {N_x}")
    print(f"- Number of distinct vertical placements (N_y) = {N_y}")
    
    print("\nPlugging these values into the formula gives:")
    print(f"D = log({N_y}) / log({s_y_inv:.2f}) + log({N_x}) / log({denominator_base})")
    print("Since log(1) = 0, the equation simplifies to:")
    
    # Final equation and result
    print("\nFinal Equation:")
    print(f"D = log({numerator}) / log({denominator_base})")
    print(f"D = {dimension}")

calculate_fractal_dimension()