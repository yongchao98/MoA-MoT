import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension for the described fractal piano.
    """
    # Number of self-similar copies
    N = 5

    # Original keyboard dimensions
    W = 3.0
    H = 1.0

    # Black key dimensions
    # White keys are equally sized (7 keys in 3 units width) -> white key width = 3/7
    # Black keys are half as wide -> black key width = (3/7)/2 = 3/14
    w_b = 3.0 / 14.0
    # Black keys cover 9/14 of the height
    h_b = 9.0 / 14.0

    # Scaling factors
    r_x = w_b / W  # (3/14) / 3 = 1/14
    r_y = h_b / H  # (9/14) / 1 = 9/14
    
    # In the dimension formula D = log(N * r_y/r_x) / log(1/r_x),
    # the numerator part inside the log is N * r_y / r_x
    numerator_val = N * (r_y / r_x)
    
    # The denominator part inside the log is 1 / r_x
    denominator_val = 1 / r_x

    # Calculate the dimension
    dimension = math.log(numerator_val) / math.log(denominator_val)

    print("The Minkowski–Bouligand dimension (D) is calculated using the formula:")
    print("D = log(N * r_y / r_x) / log(1 / r_x)")
    print("\nSubstituting the derived values:")
    # The problem asks to output each number in the final equation.
    # N=5, r_y/r_x=9, 1/r_x=14
    print(f"D = log({N} * {int(r_y/r_x)}) / log({int(1/r_x)})")
    print(f"D = log({int(numerator_val)}) / log({int(denominator_val)})")
    print(f"\nThe calculated dimension is: {dimension}")

if __name__ == "__main__":
    calculate_fractal_dimension()