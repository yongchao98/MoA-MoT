import math

def solve_modified_logistic_map():
    """
    Calculates the equilibrium point of a modified logistic map for a given R.
    The standard logistic map is X_n+1 = R * X_n * (1 - X_n).
    A proposed modification is: X_n+1 = R * X_n * (1 - X_n) + R - 2*X_n.
    
    To find the equilibrium point X*, we solve X* = f(X*, R):
    X = R*X*(1-X) + R - 2*X
    This simplifies to the quadratic equation: R*X^2 + (3-R)*X - R = 0.
    """
    R = 3.57

    # Coefficients for the quadratic equation a*X^2 + b*X + c = 0
    a = R
    b = 3 - R
    c = -R

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        print("The equation has no real equilibrium points.")
        return

    # Calculate the two possible solutions for X
    # We are interested in the positive solution
    x1 = (-b + math.sqrt(discriminant)) / (2*a)
    x2 = (-b - math.sqrt(discriminant)) / (2*a)
    
    equilibrium_point = x1 if x1 > 0 else x2

    # Print the modified equation and the result
    # The final code needs to output each number in the final equation.
    # Modified Map: Xn+1 = R * Xn * (1 - Xn) + R - 2*Xn
    # Equilibrium equation: R*X^2 + (3-R)*X - R = 0
    # For R = 3.57, the equation is 3.57*X^2 + (3-3.57)*X - 3.57 = 0
    
    print("Modified Logistic Map Equation: X_n+1 = R * X_n * (1 - X_n) + R - 2 * X_n")
    print("\nFor R = 3.57, the equilibrium equation is:")
    print(f"{a}*X^2 + ({b})*X + ({c}) = 0")
    print(f"\nWhich simplifies to:")
    print(f"3.57*X^2 - 0.57*X - 3.57 = 0")
    
    print(f"\nThe calculated equilibrium point is approximately: {equilibrium_point}")

solve_modified_logistic_map()
<<<1.0830134453781512>>>