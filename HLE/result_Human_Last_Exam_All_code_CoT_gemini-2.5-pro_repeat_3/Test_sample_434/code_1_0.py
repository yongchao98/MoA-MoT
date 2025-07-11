import math

def find_blowup_condition():
    """
    For the given system of ODEs, we determine the range of initial values y(0)
    that lead to a blow-up solution, for a given x(0) > 1.

    The system is:
    x'(t) = -3*x(t)*y(t)
    y'(t) = -y^2(t) - x(t) + 1

    A blow-up occurs if the initial condition (x(0), y(0)) lies "below" a separatrix
    in the phase plane. The separatrix is the trajectory passing through the saddle
    point (1, 0). Its equation is found to be y^2 = 2*x + 1 - 3*x^(2/3).

    The condition for blow-up is therefore y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).

    We will demonstrate this for a specific value, x(0) = 8.
    """
    x0 = 8

    # Calculate the terms for the critical value of y(0)
    term_x0_pow = x0**(2/3)
    term1 = 2 * x0
    term2 = 1
    term3 = 3 * term_x0_pow

    # Calculate the value under the square root
    radicand = term1 + term2 - term3

    # Calculate the final critical value
    y_crit = math.sqrt(radicand)

    print(f"Given the initial condition x(0) = {x0}, we want to find the values of y(0) for which the solution blows up.")
    print("The solution blows up if y(0) is less than a critical value, y_crit.")
    print("This is derived from the separatrix equation: y^2 = 2*x + 1 - 3*x^(2/3).")
    print("\nWe calculate y_crit by substituting x = x(0) into the equation for the separatrix:")
    print(f"y_crit = sqrt(2 * x(0) + 1 - 3 * x(0)^(2/3))")
    
    # Printing each number in the final equation as it's computed
    print(f"\nFor x(0) = {x0}:")
    print(f"y_crit = sqrt(2 * {x0} + {1} - 3 * ({x0})^(2/3))")
    print(f"First, we compute the power: {x0}^(2/3) = {term_x0_pow:.4f}")
    print(f"y_crit = sqrt({term1} + {term2} - 3 * {term_x0_pow:.4f})")
    print(f"y_crit = sqrt({term1} + {term2} - {term3:.4f})")
    print(f"y_crit = sqrt({radicand:.4f})")
    print(f"y_crit = {y_crit:.4f}")

    print(f"\nTherefore, for x(0) = {x0}, the solution blows up if y(0) < {y_crit:.4f}.")
    print(f"The set of values for y(0) is (-infinity, {y_crit:.4f}), which is approximately (-inf, sqrt(5)).")
    
    # The final answer is the critical value of y(0) for x(0) = 8.
    print(f"<<<{y_crit}>>>")

find_blowup_condition()