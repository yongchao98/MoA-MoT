import math

def find_blowup_condition(x0: float):
    """
    For the given system of differential equations, this function calculates
    the condition on the initial value y(0) that leads to a finite-time
    blow-up, assuming x(0) = x0 > 1.

    Args:
        x0: The initial value for x(t), which must be greater than 1.
    """
    if x0 <= 1:
        print("Error: The analysis is valid only for x(0) > 1.")
        return

    # The condition for blow-up is y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).
    # We will calculate the terms of this expression.
    
    term_2x0 = 2 * x0
    term_1 = 1
    term_3x_pow_2_3 = 3 * (x0**(2/3))

    expression_inside_sqrt = term_2x0 + term_1 - term_3x_pow_2_3

    # For x0 > 1, the expression inside the square root is always positive.
    y_critical = math.sqrt(expression_inside_sqrt)

    print("For the system of differential equations:")
    print("  x'(t) = -3*x(t)*y(t)")
    print("  y'(t) = -y(t)^2 - x(t) + 1")
    print(f"\nWith the initial condition x(0) = {x0} (which is > 1), the solution blows up if y(0) is less than a critical value.")
    print("\nThe blow-up condition is given by the inequality:")
    print("  y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3))")

    print("\nSubstituting the value of x(0):")
    # Outputting each number in the final equation
    print(f"  y(0) < sqrt(2 * {x0} + 1 - 3 * {x0}^(2/3))")
    print(f"  y(0) < sqrt({term_2x0} + {term_1} - {term_3x_pow_2_3:.4f})")
    print(f"  y(0) < sqrt({expression_inside_sqrt:.4f})")
    print(f"  y(0) < {y_critical:.4f}")

# You can change the value of x(0) here. Let's use x(0) = 4 as an example.
initial_x = 4.0
find_blowup_condition(initial_x)
