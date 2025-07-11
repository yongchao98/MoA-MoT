import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski-Bouligand dimension for the described fractal piano.
    """
    # Parameters from the problem description
    N = 5  # Number of black keys in one octave
    
    # Original keyboard dimensions
    W_orig = 3
    H_orig = 1
    
    # Black key dimensions
    W_black = 3 / 14
    H_black = 9 / 14

    # Scaling factors
    s_x = W_black / W_orig # 1/14
    s_y = H_black / H_orig # 9/14
    
    # The dimension D is calculated using the formula:
    # D = 1 + log(N * s_y) / log(1/s_x)  --- derived from N * s_y * s_x^(D-1) = 1
    # which simplifies to D = 1 + log(45/14) / log(14)
    
    val_num_arg = 45 / 14
    val_den_arg = 14

    log_numerator = math.log(val_num_arg)
    log_denominator = math.log(val_den_arg)
    
    # D-1 term
    d_minus_1 = log_numerator / log_denominator
    
    # Final dimension D
    D = 1 + d_minus_1

    print("The Minkowski–Bouligand dimension D is the solution to the equation:")
    print("5 * (9/14) * (1/14)^(D-1) = 1")
    print("\nSolving for D gives the formula:")
    print(f"D = 1 + log({45}/{14}) / log({14})")
    print("\nCalculating the values:")
    print(f"D = 1 + log({val_num_arg:.4f}) / log({val_den_arg})")
    print(f"D = 1 + {log_numerator:.4f} / {log_denominator:.4f}")
    print(f"D = 1 + {d_minus_1:.4f}")
    print(f"\nThe final dimension is D = {D:.4f}")

calculate_fractal_dimension()
<<<1.4424>>>