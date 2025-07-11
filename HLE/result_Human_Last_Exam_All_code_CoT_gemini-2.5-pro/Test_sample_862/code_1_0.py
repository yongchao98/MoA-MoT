import numpy as np
from scipy.optimize import root_scalar

def solve_for_C():
    """
    Solves for the constant C based on the derived transcendental equation.
    """
    # The transcendental equation for theta is tan(2*theta) = 2*(3*pi/4 - theta).
    # We define a function whose root we need to find.
    # The variable 'x' in the function corresponds to theta in the equation.
    equation = lambda x: np.tan(2 * x) - 2 * (3 * np.pi / 4 - x)

    # We search for a root in the interval (pi/4, pi/2), as derived in the logic.
    pi_val = np.pi
    # Add a small epsilon to avoid the boundaries where tan is undefined.
    bracket = (pi_val / 4 + 1e-9, pi_val / 2 - 1e-9)
    
    try:
        sol = root_scalar(equation, bracket=bracket)
        if sol.converged:
            theta = sol.root
            C = np.tan(theta)**2
            
            print(f"The problem reduces to solving a transcendental equation for an angle theta.")
            print(f"The equation is: tan(2 * theta) = 2 * (3*pi/4 - theta)")
            print(f"The equation can be written with all terms on one side:")
            print(f"tan(2 * {theta:.5f}) - 2 * (3*{pi_val:.5f}/4 - {theta:.5f}) = {equation(theta):.5f}")
            print("\nSolving this equation numerically for theta in the interval (pi/4, pi/2):")
            print(f"theta = {theta:.8f} radians")
            
            print("\nThe constant C is given by tan(theta)^2.")
            print(f"C = tan({theta:.8f})^2")
            print(f"C = {np.tan(theta):.8f}^2")
            print(f"C = {C:.8f}")

        else:
            print("Numerical solver did not converge.")
    except ValueError as e:
        print(f"An error occurred during root finding: {e}")
        print("This might happen if the function values at the bracket endpoints have the same sign.")

solve_for_C()
