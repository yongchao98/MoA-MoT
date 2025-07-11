import math

def solve_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension for the described fractal piano.
    """
    # Parameters from the problem description
    N = 5  # Number of black keys in an octave
    original_width = 3
    original_height = 1
    
    # Dimensions of a single black key
    black_key_width = (1/2) * (original_width / 7)
    black_key_height = (9/14) * original_height
    
    # Anisotropic scaling factors
    r_x = black_key_width / original_width
    r_y = black_key_height / original_height
    
    print("Step 1: Determine the fractal parameters.")
    print(f"The fractal is composed of N = {N} smaller copies.")
    print(f"The horizontal scaling factor, r_x = ({black_key_width:.4f}) / {original_width} = {r_x:.4f}")
    print(f"The vertical scaling factor, r_y = ({black_key_height:.4f}) / {original_height} = {r_y:.4f}")
    print("-" * 20)
    
    print("Step 2: Apply the formula for the box-counting dimension (D).")
    print("The formula for this type of self-affine set is: D = log(N * r_y / r_x) / log(1/r_x)")
    print("-" * 20)

    print("Step 3: Calculate the result.")
    # Arguments for the log functions
    numerator_arg = N * r_y / r_x
    denominator_arg = 1 / r_x
    
    # Calculate the dimension
    dimension = math.log(numerator_arg) / math.log(denominator_arg)
    
    print("Substituting the values into the formula:")
    print(f"D = log({N} * {r_y:.4f} / {r_x:.4f}) / log(1 / {r_x:.4f})")
    print("This simplifies to the final equation:")
    print(f"D = log({int(numerator_arg)}) / log({int(denominator_arg)})")
    print(f"D = {math.log(numerator_arg):.4f} / {math.log(denominator_arg):.4f}")
    print(f"\nThe calculated dimension is D = {dimension:.4f}")

solve_fractal_dimension()