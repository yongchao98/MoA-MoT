import math

def calculate_piano_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the described fractal piano.
    """
    # Step 1: Define the parameters based on the problem description.

    # N is the number of self-similar copies at each iteration.
    # This corresponds to the number of black keys in one octave.
    N = 5

    # s_x is the scaling factor in the x-direction (width).
    # The keyboard is 3 units wide with 7 equally sized white keys.
    # Width of a white key = 3 / 7
    # Width of a black key = 1/2 * (width of a white key) = 1/2 * (3/7) = 3/14
    # The scaling factor s_x is the ratio of the new width to the old width.
    # s_x = (width of black key) / (width of keyboard) = (3/14) / 3
    s_x = 1/14

    # s_y is the scaling factor in the y-direction (height).
    # The keyboard is 1 unit high.
    # Height of a black key = 9/14 of the total height = 9/14
    # The scaling factor s_y is the ratio of the new height to the old height.
    # s_y = (height of black key) / (height of keyboard) = (9/14) / 1
    s_y = 9/14

    # Step 2: Use the formula for the box-counting dimension of a self-affine set.
    # The general formula is D = (log(N) + log(s_y) - log(s_x)) / log(1/s_x).
    # Let's simplify this for our values:
    # D = log(N * s_y / s_x) / log(1/s_x)
    # D = log(5 * (9/14) / (1/14)) / log(1 / (1/14))
    # D = log(5 * 9) / log(14)
    # D = log(45) / log(14)

    numerator_val = 45
    denominator_val = 14

    # Step 3: Calculate the numerical value of the dimension.
    log_numerator = math.log(numerator_val)
    log_denominator = math.log(denominator_val)
    dimension = log_numerator / log_denominator

    # Step 4: Print the explanation and the final equation with all its components.
    print("The Minkowski–Bouligand dimension (D) is calculated as follows:")
    print("\n1. Parameters:")
    print(f"   - Number of copies (N): {N}")
    print(f"   - Scaling factor for width (s_x): 1/14")
    print(f"   - Scaling factor for height (s_y): 9/14")

    print("\n2. Formula:")
    print("   The dimension D for this self-affine set simplifies to:")
    print(f"   D = log({numerator_val}) / log({denominator_val})")

    print("\n3. Calculation:")
    print(f"   log({numerator_val}) = {log_numerator}")
    print(f"   log({denominator_val}) = {log_denominator}")
    print(f"   D = {log_numerator} / {log_denominator}")

    print("\n4. Final Answer:")
    print(f"   The dimension D is: {dimension}")

calculate_piano_fractal_dimension()