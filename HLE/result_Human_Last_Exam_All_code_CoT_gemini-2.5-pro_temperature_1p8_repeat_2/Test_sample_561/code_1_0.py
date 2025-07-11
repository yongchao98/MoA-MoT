import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension of the fractal piano keyboard.
    """
    # Define parameters from the problem description
    keyboard_width = 3.0
    keyboard_height = 1.0
    num_white_keys = 7
    # N is the number of self-similar copies, which is the number of black keys
    N = 5

    # Calculate the dimensions of a single black key
    white_key_width = keyboard_width / num_white_keys
    black_key_width = white_key_width / 2.0
    black_key_height = keyboard_height * (9.0 / 14.0)

    # Calculate the scaling factors in x and y directions
    s_x = black_key_width / keyboard_width
    s_y = black_key_height / keyboard_height
    
    # For a self-affine fractal of this type, the dimension D is given by:
    # D = log(N * s_y / s_x) / log(1 / s_x)
    # This formula simplifies when we substitute the values.
    # N = 5
    # s_y / s_x = (9/14) / (1/14) = 9
    # 1 / s_x = 1 / (1/14) = 14
    # So, D = log(5 * 9) / log(14) = log(45) / log(14)

    numerator_simplified = 45
    denominator_simplified = 14
    
    # Calculate the final dimension
    dimension = math.log(numerator_simplified) / math.log(denominator_simplified)

    # Print the explanation and the final equation with its components
    print("The Minkowski-Bouligand dimension (D) for this fractal is calculated from its self-affine properties.")
    print("\nNumber of copies (N): 5")
    print("Scaling factor in x (s_x): 1/14")
    print("Scaling factor in y (s_y): 9/14")
    
    print("\nThe dimension formula simplifies to D = log(N * (s_y/s_x)) / log(1/s_x).")
    print("This gives the final equation:")
    print(f"D = log({numerator_simplified}) / log({denominator_simplified})")
    
    print("\nThe calculated dimension is:")
    print(dimension)


if __name__ == '__main__':
    calculate_fractal_dimension()