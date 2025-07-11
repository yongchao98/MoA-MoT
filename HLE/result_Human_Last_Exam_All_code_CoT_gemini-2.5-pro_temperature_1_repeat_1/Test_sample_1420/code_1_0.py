import math

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α)
    using the epsilon expansion for d=3.
    """
    # Step 1: Define the parameters for the calculation.
    # d is the spatial dimension.
    d = 3
    # d_c is the upper critical dimension for an Ising-like system.
    d_c = 4
    # n is the number of components of the order parameter.
    # For the Ising model, n=1.
    n = 1

    print("Calculating the specific heat scaling exponent alpha (α) using the epsilon expansion.")
    print(f"The calculation is for a system with spatial dimension d = {d}.")
    print(f"The upper critical dimension for this universality class is d_c = {d_c}.")
    print(f"We assume the Ising model, so the order parameter has n = {n} component(s).\n")

    # Step 2: Calculate epsilon (ε).
    epsilon = d_c - d
    print(f"First, we calculate the expansion parameter ε (epsilon):")
    print(f"ε = d_c - d = {d_c} - {d} = {epsilon}\n")

    # Step 3: Use the one-loop formula for alpha and substitute values.
    print("Next, we use the one-loop renormalization group (RG) formula for α:")
    print(f"α = (4 - n) / (2 * (n + 8)) * ε\n")
    
    print("Now, we substitute the values n=1 and ε=1 into the equation:")
    
    # Perform the calculation, showing each step.
    numerator = 4 - n
    intermediate_denominator = n + 8
    denominator = 2 * intermediate_denominator
    
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / (2 * {intermediate_denominator}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    
    alpha = (numerator / denominator) * epsilon
    
    print(f"α = {numerator/denominator} * {epsilon}")
    print(f"α = {alpha}\n")

    print(f"The result to first order in the epsilon expansion is α = 1/6, which is approximately {alpha:.3f}.")

calculate_alpha_exponent()