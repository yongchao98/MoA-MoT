import math

def calculate_blowup_condition(x0):
    """
    Calculates the critical value for y(0) that determines blow-up behavior
    for a given x(0) > 1.
    """
    # Check if the term inside the square root is non-negative.
    # The function g(x) = 2x + 1 - 3x^(2/3) is positive for x > 1.
    g_x0 = 2*x0 + 1 - 3*x0**(2/3)
    if g_x0 < 0:
        print(f"Error: The value x(0)={x0} is not valid as it leads to a complex threshold.")
        return

    # Calculate the components of the expression
    term1 = 2 * x0
    term2 = 1
    term3 = 3 * (x0**(2/3))

    # Calculate the value inside the square root
    y_crit_sq = term1 + term2 - term3
    
    # Calculate the critical threshold
    y_crit = math.sqrt(y_crit_sq)

    # Print the explanation and result
    print("The system blows up if y(0) is less than a critical value derived from the separatrix equation.")
    print(f"The separatrix is y^2 = 2x + 1 - 3x^(2/3).")
    print(f"For an initial condition x(0) = {x0}, the critical value for y(0) is:")
    print(f"y_crit = sqrt(2*x(0) + 1 - 3*x(0)^(2/3))")
    print(f"         = sqrt(2*{x0} + {term2} - 3*({x0})^(2/3))")
    print(f"         = sqrt({term1} + {term2} - {term3:.4f})")
    print(f"         = sqrt({y_crit_sq:.4f})")
    print(f"         = {y_crit:.4f}")
    print("\nTherefore, the solution blows up if y(0) is less than this value.")
    print(f"Blow-up occurs for y(0) < {y_crit:.4f}")

# You can change this value to any number greater than 1
x_initial = 8.0
calculate_blowup_condition(x_initial)
