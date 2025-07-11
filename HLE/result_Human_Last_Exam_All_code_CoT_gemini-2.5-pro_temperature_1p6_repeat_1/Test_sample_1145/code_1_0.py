import math

def calculate_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a given interval.
    """
    # Given parameters
    n = 2
    x1_div_a = 0.495
    x2_div_a = 0.505

    # The probability is found by integrating the probability density |ψ(x)|² from x1 to x2.
    # The analytical solution for P = ∫[x1 to x2] (2/a)sin²(nπx/a) dx is:
    # P = [x/a - (1/(2nπ))sin(2nπx/a)] from x1 to x2.
    # P = (x2/a - x1/a) - (1/(2nπ)) * [sin(2nπx2/a) - sin(2nπx1/a)]

    # Calculate the components of the equation
    delta_x_div_a = x2_div_a - x1_div_a
    
    coeff_inv = 2 * n * math.pi
    coeff = 1 / coeff_inv

    arg2 = 2 * n * math.pi * x2_div_a
    arg1 = 2 * n * math.pi * x1_div_a
    
    sin_term2 = math.sin(arg2)
    sin_term1 = math.sin(arg1)
    
    # Calculate the final probability
    probability = delta_x_div_a - coeff * (sin_term2 - sin_term1)

    # Print the equation with numerical values as requested
    print("The probability P is given by the formula:")
    print("P = (x₂/a - x₁/a) - (1/(2nπ)) * [sin(2nπx₂/a) - sin(2nπx₁/a)]")
    print("\nWith n = 2, x₁/a = 0.495, and x₂/a = 0.505, the equation becomes:")
    print(f"P = ({x2_div_a} - {x1_div_a}) - (1 / (2 * {n} * π)) * [sin(2 * {n} * π * {x2_div_a}) - sin(2 * {n} * π * {x1_div_a})]")
    print(f"P = {delta_x_div_a} - (1 / {coeff_inv:.6f}) * [{sin_term2:.6f} - ({sin_term1:.6f})]")
    print(f"P = {delta_x_div_a} - {coeff:.6f} * [{sin_term2 - sin_term1:.6f}]")

    # Print the final result
    print("\nThe final calculated probability is:")
    print(f"P = {probability:.8f}")

    return probability

# Run the calculation and store the result
final_probability = calculate_probability()

# The final answer in the requested format will be outputted at the end.
# Note: The printing is done inside the function.
