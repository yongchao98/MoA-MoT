import math

def find_blowup_condition(x0: float):
    """
    For the given system of differential equations, this function determines the
    range of initial values y(0) that lead to a blow-up solution, given an
    initial value x(0) > 1.

    Args:
        x0: The initial value x(0), which must be greater than 1.
    """
    if x0 <= 1:
        print("Error: This analysis is valid for x(0) > 1.")
        return

    # The condition for blow-up is determined by a separatrix curve, which is the
    # stable manifold of the saddle point at (1, 0).
    # The equation for this separatrix is y^2 = c1*x + c2 - c3*x^(p1/p2).
    c1 = 2
    c2 = 1
    c3 = 3
    p1 = 2
    p2 = 3

    # Calculate the term under the square root for the critical y value
    radicand = c1 * x0 + c2 - c3 * x0**(p1 / p2)
    
    # For x0 > 1, the radicand is always non-negative.
    y_critical = math.sqrt(radicand)

    print(f"For the initial condition x(0) = {x0}, we determine the values of y(0) for which the solution blows up.")
    print("\nThe analysis shows that a blow-up occurs if the initial point (x(0), y(0)) lies below a specific curve (a separatrix).")
    print("This separatrix is the stable manifold of the saddle point (1,0) and has the equation:")
    print(f"y^2 = {c1}*x + {c2} - {c3}*x^({p1}/{p2})\n")
    
    print("To find the critical value for y(0), we substitute x = x(0) into this equation:")
    print(f"y_crit^2 = {c1}*({x0}) + {c2} - {c3}*({x0})^({p1}/{p2})")
    print(f"y_crit^2 = {c1*x0} + {c2} - {c3*x0**(p1/p2):.4f}")
    print(f"y_crit^2 = {radicand:.4f}")
    print(f"y_crit = {y_critical:.4f}\n")
    
    print("The solution blows up if y(0) is less than this critical value.")
    print(f"Therefore, the condition for blow-up is: y(0) < {y_critical:.4f}")

# Set a value for x(0) > 1. For example, x(0) = 2.
initial_x = 2.0
find_blowup_condition(initial_x)