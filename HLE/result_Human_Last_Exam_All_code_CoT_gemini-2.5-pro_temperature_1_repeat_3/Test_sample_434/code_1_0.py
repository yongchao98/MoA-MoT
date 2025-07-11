import math

def print_blowup_condition(x0):
    """
    For a given initial condition x(0) = x0 > 1, this function prints the
    condition on y(0) for which the solution to the system of ODEs blows up.

    The system is:
    x'(t) = -3*x(t)*y(t)
    y'(t) = -y(t)**2 - x(t) + 1
    """
    if x0 <= 1:
        print("Error: This analysis is valid only for x(0) > 1.")
        return

    # The critical value for y(0) is derived from the separatrix equation:
    # y^2 = 2*x + 1 - 3*x**(2/3)
    # The solution blows up if the initial condition (x0, y0) is "below"
    # the upper branch of this separatrix curve.
    
    # Coefficients of the equation y^2 = c2*x + c1 + c3*x^exp
    c2 = 2
    c1 = 1
    c3 = -3
    exp = 2/3

    # For x0 > 1, the term 2*x0 + 1 - 3*x0**(2/3) is always positive,
    # so its square root is a real number.
    rhs_squared = c2 * x0 + c1 + c3 * (x0**exp)
    y_crit = math.sqrt(rhs_squared)

    print(f"For the given initial condition x(0) = {x0}:")
    print("The solution to the system blows up if y(0) is less than a critical value.")
    print("\nThis critical value is derived from the separatrix equation:")
    # Using f-string formatting to show the equation structure
    print(f"y(0)^2 = ({c2})*x(0) + ({c1}) + ({c3})*x(0)^({exp:.2f})")
    
    print(f"\nTo find the critical value, we evaluate the condition y(0)^2 < 2*x(0) + 1 - 3*x(0)^(2/3):")
    # Showing each number/term in the calculation
    term1 = c2 * x0
    term2 = c1
    term3_val = c3 * x0**exp
    print(f"y(0)^2 < {term1} + {term2} + ({term3_val:.4f})")
    print(f"y(0)^2 < {rhs_squared:.4f}")
    
    print("\nThe final condition for blow-up is therefore:")
    print(f"y(0) < {y_crit:.4f}")

# You can test this with any initial condition x(0) > 1.
# Let's use x(0) = 4 as an example.
x_initial = 4
print_blowup_condition(x_initial)