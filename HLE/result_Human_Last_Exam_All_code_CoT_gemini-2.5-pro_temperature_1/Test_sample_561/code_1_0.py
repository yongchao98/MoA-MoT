import math

def calculate_fractal_dimension():
    """
    Calculates and explains the Minkowski-Bouligand dimension of the fractal piano keys.
    """
    
    # Step 1: Define the parameters of the fractal construction.
    # N is the number of self-similar copies. There are 5 black keys in an octave.
    N = 5
    
    # s_x and s_y are the scaling factors in the x and y directions.
    # A new keyboard (3x1) is scaled to fit onto a black key (3/14 x 9/14).
    s_x = (3/14) / 3  # = 1/14
    s_y = (9/14) / 1  # = 9/14
    
    print("This is a self-affine fractal. We first determine its parameters:")
    print(f"Number of copies (N): {N}")
    print(f"Scaling factor in x-direction (s_x): 1/14")
    print(f"Scaling factor in y-direction (s_y): 9/14")
    print("-" * 20)

    # Step 2: State and solve the equation for the dimension D.
    # For a self-affine fractal of this type, the dimension D is the solution to:
    # N * s_y * s_x^(D-1) = 1
    # 5 * (9/14) * (1/14)^(D-1) = 1
    # Solving for D gives: D = log(N * s_y * s_x^(-1)) / log(s_x^(-1))
    # which simplifies to D = log(N * s_y / s_x) / log(1/s_x)
    # D = log(5 * (9/14) / (1/14)) / log(14)
    # D = log(5 * 9) / log(14)
    # D = log(45) / log(14)
    
    numerator = N * s_y / s_x
    denominator = 1/s_x
    
    final_dimension = math.log(numerator) / math.log(denominator)

    print("The dimension D is found by solving the equation: N * s_y * (s_x)^(D-1) = 1")
    print("Plugging in the values, this simplifies to the final equation:")
    print(f"D = log({int(numerator)}) / log({int(denominator)})")
    print("-" * 20)
    print("The calculated Minkowski-Bouligand dimension is:")
    print(final_dimension)

calculate_fractal_dimension()
<<<1.4424349448883344>>>