import math

def print_blowup_condition():
    """
    This function prints the derived condition for the solution to blow up.
    """
    # The separatrix curve is found to be y^2 = 2*x + 1 - 3*x^(2/3).
    # Blow-up occurs for initial conditions (x(0), y(0)) that lie below this curve.
    
    C1 = 2
    C2 = 1
    C3 = -3
    P1 = 2
    P2 = 3
    
    print("For the system of differential equations with the initial condition x(0) > 1,")
    print("the solution (x(t), y(t)) blows up in finite time if and only if y(0) satisfies the following inequality:")
    
    # Printing the final equation with all its numbers
    print(f"\ny(0) < sqrt({C1}*x(0) + {C2} + {C3}*x(0)**({P1}/{P2}))\n")

    # The expression under the square root is g(x) = 2x + 1 - 3x^(2/3).
    # g(1) = 0 and g'(x) = 2 - 2x^(-1/3) > 0 for x > 1.
    # So the right hand side is always a positive real number for x(0) > 1.
    print("This condition correctly implies that any non-positive y(0) (i.e., y(0) <= 0) will lead to a blow-up.")

    # We can also demonstrate this with a numerical example.
    x0_example = 4.0
    # For x=4, x^(2/3) is approx 2.5198
    y_critical = math.sqrt(C1 * x0_example + C2 + C3 * (x0_example**(P1/P2)))
    print(f"\nFor instance, if x(0) = {int(x0_example)}, the solution blows up for y(0) < {y_critical:.4f}.")

print_blowup_condition()