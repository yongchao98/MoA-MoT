import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the fractal piano keys.
    """
    # Step 1 & 2: Define the parameters of the IFS from the problem description.
    # Number of black keys in one octave. This is N, the number of transformations.
    N = 5
    
    # Dimensions of the main keyboard
    keyboard_width = 3.0
    keyboard_height = 1.0
    
    # Number of white keys
    num_white_keys = 7
    
    # Properties of the black keys
    black_key_height_fraction = 9.0 / 14.0
    black_key_width_fraction = 0.5  # half as wide as white keys

    # Step 3: Calculate the scaling factors s_x and s_y.
    # First, find the dimensions of a single black key.
    white_key_width = keyboard_width / num_white_keys
    black_key_width = white_key_width * black_key_width_fraction
    black_key_height = keyboard_height * black_key_height_fraction

    # The transformation scales the whole keyboard to fit onto a black key.
    s_x = black_key_width / keyboard_width
    s_y = black_key_height / keyboard_height

    print(f"Number of recursive components (N): {N}")
    print(f"Horizontal scaling factor (s_x): {s_x:.4f}")
    print(f"Vertical scaling factor (s_y): {s_y:.4f}")
    print("-" * 20)

    # Step 4: Apply the box-counting dimension formula for this type of self-affine set.
    # The formula is D = log(N * s_y / s_x) / log(1 / s_x)
    
    # Calculate the terms in the formula.
    ratio_of_scalings = s_y / s_x
    numerator_term = N * ratio_of_scalings
    denominator_term = 1 / s_x

    print("The dimension D is calculated using the formula:")
    print("D = log(N * s_y / s_x) / log(1 / s_x)\n")
    
    print("Let's calculate the values for the final equation:")
    print(f"The term in the numerator is N * s_y / s_x = {N} * ({s_y:.4f} / {s_x:.4f}) = {numerator_term:.1f}")
    print(f"The term in the denominator is 1 / s_x = 1 / {s_x:.4f} = {denominator_term:.1f}")
    
    print("\nThe final equation is:")
    print(f"D = log({numerator_term:.1f}) / log({denominator_term:.1f})")

    # Step 5: Compute the final value.
    dimension = math.log(numerator_term) / math.log(denominator_term)
    
    print(f"\nThe calculated Minkowski–Bouligand dimension is: {dimension}")

calculate_fractal_dimension()