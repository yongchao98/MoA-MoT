import math

def calculate_piano_fractal_dimension():
    """
    This function calculates the Minkowski-Bouligand dimension of the fractal
    described in the problem.
    """

    # Step 1: Define the geometry of the keyboard to find scaling factors.
    
    # An octave keyboard has 7 white keys and 5 black keys.
    # The total width of the keyboard is 3 units, covered by the 7 white keys.
    keyboard_width = 3.0
    # The width of a single white key.
    white_key_width = keyboard_width / 7.0
    # A black key is half as wide as a white key.
    black_key_width = white_key_width / 2.0

    # The total height of the keyboard is 1 unit.
    keyboard_height = 1.0
    # The black keys cover 9/14 of the keyboard's height.
    black_key_height = keyboard_height * (9.0 / 14.0)

    # Step 2: Determine the parameters for the fractal dimension formula.
    
    # N: The number of self-similar copies. A new keyboard is placed on each
    # of the 5 black keys, so N = 5.
    N = 5

    # r_x: The horizontal scaling factor. The whole keyboard (width=3) is mapped
    # onto a single black key (width=3/14).
    # r_x = (3/14) / 3 = 1/14
    r_x = black_key_width / keyboard_width

    # r_y: The vertical scaling factor. The whole keyboard (height=1) is mapped
    # onto a single black key (height=9/14).
    # r_y = (9/14) / 1 = 9/14
    r_y = black_key_height / keyboard_height

    # Since r_x (1/14) is smaller than r_y (9/14), the horizontal direction
    # is the "faster" shrinking direction.
    
    # Step 3: Formulate the equation for the dimension D.
    # The formula is D = log(N * r_y / r_x) / log(1 / r_x)
    
    # Let's calculate the integer values for the arguments of the logarithms.
    numerator_arg = N * (r_y / r_x)
    denominator_arg = 1 / r_x
    
    # We can represent these values exactly with integers.
    # numerator = 5 * ( (9/14) / (1/14) ) = 5 * 9 = 45
    # denominator = 1 / (1/14) = 14
    final_numerator_val = 45
    final_denominator_val = 14
    
    # Step 4: Calculate the final dimension.
    dimension = math.log(final_numerator_val) / math.log(final_denominator_val)

    # Print the final result, including the equation.
    print(f"The fractal dimension D is given by the formula D = log(A) / log(B)")
    print(f"Where A = N * (vertical scaling / horizontal scaling) = 5 * ( (9/14) / (1/14) ) = {final_numerator_val}")
    print(f"And B = 1 / (horizontal scaling) = 1 / (1/14) = {final_denominator_val}")
    print("\n--- Final Equation ---")
    print(f"D = log({final_numerator_val}) / log({final_denominator_val})")
    print(f"\nThe Minkowski-Bouligand dimension is: {dimension}")

calculate_piano_fractal_dimension()