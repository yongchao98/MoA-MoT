import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension of the fractal piano.
    """
    # N is the number of self-similar copies at each iteration.
    # In one octave, there are 5 black keys, so N=5.
    N = 5

    # r is the scaling factor.
    # The main keyboard (3x1) must be scaled to fit on a black key.
    # Width of a black key = (3/7)/2 = 3/14
    # Height of a black key = 1 * (9/14) = 9/14
    # To maintain aspect ratio, we scale the main keyboard by 'r'.
    # Scaled width: r * 3 <= 3/14  => r <= 1/14
    # Scaled height: r * 1 <= 9/14
    # The more restrictive factor is r = 1/14.
    r = 1.0 / 14.0

    # For the dimension formula D = log(N) / log(1/r), we need 1/r.
    reciprocal_r = 1.0 / r

    # Calculate the dimension.
    dimension = math.log(N) / math.log(reciprocal_r)

    # Print the explanation and the result.
    print("To find the Minkowski–Bouligand dimension (D), we use the formula:")
    print("D = log(N) / log(1/r)")
    print("\nWhere:")
    print(f"- N is the number of self-similar copies per iteration, which is the number of black keys: {N}")
    print(f"- r is the scaling factor. We determined r = 1/14, so 1/r = {int(reciprocal_r)}")
    print("\nThe final equation is:")
    print(f"D = log({N}) / log({int(reciprocal_r)})")
    print(f"\nThe calculated dimension is: {dimension}")
    print(f"\n<<<{dimension}>>>")

calculate_fractal_dimension()