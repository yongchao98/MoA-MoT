import math

def get_blowup_y0_bound(x0):
    """
    For the given system of differential equations, this function calculates the
    upper bound for the initial condition y(0) that leads to a blow-up solution,
    assuming the initial condition x(0) = x0 is greater than 1.

    Args:
        x0 (float): The initial value x(0), which must be > 1.

    Returns:
        float: The upper bound for y(0). Blow-up occurs if y(0) is less than this value.
    
    Raises:
        ValueError: If x0 <= 1, as the analysis is valid only for x0 > 1.
    """
    if x0 <= 1:
        # The derivation relies on x0 > 1, which ensures the term under the square root is positive.
        raise ValueError("This analysis assumes x(0) > 1.")
    
    # These are the coefficients from the derived separatrix equation: y^2 = 2*x + 1 - 3*x^(2/3)
    coeff_x = 2
    constant = 1
    coeff_x_pow = -3
    power_num = 2
    power_den = 3
    
    # Calculate the expression under the square root
    value_under_sqrt = coeff_x * x0 + constant + coeff_x_pow * (x0**(power_num / power_den))
    
    # The blow-up condition is y(0) < sqrt(value_under_sqrt)
    y0_upper_bound = math.sqrt(value_under_sqrt)
    
    return y0_upper_bound

# --- Main execution ---
# The question asks for the values of y(0) for which the solution blows up.
# This depends on the value of x(0), which is given as x(0) > 1.
# The blow-up occurs if the initial point (x(0), y(0)) is "below" a critical curve (a separatrix).

# The equation for the separatrix was found to be y^2 = 2*x + 1 - 3*x^(2/3).
# Blow-up occurs for any y(0) <= 0, and for y(0) > 0 if the point is below this curve.
# This gives the combined condition: y(0) < sqrt(2*x(0) + 1 - 3*x(0)**(2/3)).

print("Based on the analysis of the system's phase portrait, a blow-up occurs if the initial")
print("condition y(0) is below a critical value determined by x(0).")
print("\nThe final equation for the blow-up condition is derived from the separatrix y^2 = 2*x + 1 - 3*x^(2/3).")
print("Each number in the final equation is presented below:")

final_equation_str = "y(0) < sqrt({0}*x(0) + {1} - {2}*x(0)**({3}/{4}))"
print(final_equation_str.format(2, 1, 3, 2, 3))

print("\n--- Example Calculation ---")
# The user can modify this example value for x(0).
x0_example = 8.0

try:
    # Calculate the bound for the example x0
    y0_upper_bound = get_blowup_y0_bound(x0_example)
    
    print(f"For a given initial condition x(0) = {x0_example}, the solution blows up if:")
    print(f"y(0) < {y0_upper_bound:.4f}")
    print(f"Thus, the range for y(0) is (-infinity, {y0_upper_bound:.4f}).")
    
except ValueError as e:
    print(e)