import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension of the fractal piano keys.
    """
    # Step 1: Define the problem parameters from the description.
    # N is the number of self-similar copies at each iteration.
    N = 5  # There are 5 black keys in an octave.

    # Initial keyboard dimensions.
    W0 = 3.0
    H0 = 1.0

    # Properties of the keys.
    num_white_keys = 7
    # A black key is half as wide as a white key.
    black_key_width_ratio = 0.5
    # A black key covers 9/14 of the keyboard height.
    black_key_height_num = 9.0
    black_key_height_den = 14.0

    # Step 2: Calculate the scaling factors for the affine transformation.
    # First, find the dimensions of a single black key.
    width_white_key = W0 / num_white_keys
    width_black_key = width_white_key * black_key_width_ratio
    height_black_key = H0 * (black_key_height_num / black_key_height_den)

    # The scaling factors are the ratio of the new dimensions (black key)
    # to the old dimensions (full keyboard).
    s_x = width_black_key / W0
    s_y = height_black_key / H0

    # Step 3: Identify the smaller and larger scaling factors.
    s_min = min(s_x, s_y)
    s_max = max(s_x, s_y)

    # Step 4: Calculate the terms for the dimension formula.
    # The dimension D is given by: D = log(N * (s_max / s_min)) / log(1 / s_min)
    # Let's calculate the arguments of the logarithms.
    numerator_arg = N * (s_max / s_min)
    denominator_arg = 1 / s_min

    # Step 5: Compute the final dimension.
    dimension = math.log(numerator_arg) / math.log(denominator_arg)

    # Step 6: Print the final equation and the result.
    print("The Minkowski-Bouligand dimension (D) for this self-affine fractal is calculated as follows:")
    print(f"D = log(N * (s_max / s_min)) / log(1 / s_min)")
    print(f"D = log({N} * (({s_y:.4f}) / ({s_x:.4f}))) / log(1 / ({s_x:.4f}))")
    print(f"D = log({N} * {s_max/s_min:.0f}) / log({1/s_min:.0f})")
    print("\nThe final equation for the dimension is:")
    print(f"D = log({int(numerator_arg)}) / log({int(denominator_arg)})")
    print(f"\nThe calculated dimension is: {dimension:.4f}")

calculate_fractal_dimension()