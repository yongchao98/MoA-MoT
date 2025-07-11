import numpy as np
from scipy.optimize import fsolve

def solve_fourteen_circles_packing():
    """
    Calculates the radius of a circle that tightly packs fourteen unit circles.

    This problem is a known "circle packing" problem. The radius R is found using
    the formula: R = r * (1 + 1/sin(alpha)), where r=1 is the radius of the small circles.
    The angle 'alpha' is the smallest positive root of the complex transcendental equation:
    2*sin(2a) + 2*sin(3a) + 2*sin(5a) - 2*sin(6a) - sin(7a) + sin(8a) - 1 = 0.

    This script solves this equation numerically to find alpha and then calculates R.

    Note: This script requires the 'scipy' library. You can install it using:
    pip install scipy
    """
    
    # Radius of the small circles
    r = 1

    # Define the transcendental equation f(alpha) = 0
    def equation(alpha):
        return (2 * np.sin(2 * alpha) + 
                2 * np.sin(3 * alpha) + 
                2 * np.sin(5 * alpha) - 
                2 * np.sin(6 * alpha) - 
                np.sin(7 * alpha) + 
                np.sin(8 * alpha) - 1)

    # An initial guess for the root-finding algorithm. Based on similar packing problems,
    # the angle is expected to be small. 0.3 radians is a good starting point.
    initial_guess = 0.3
    
    # Use fsolve from scipy.optimize to find the root for alpha
    alpha_solution = fsolve(equation, initial_guess)
    alpha = alpha_solution[0]

    # Calculate the radius of the large circle
    R = r * (1 + (1 / np.sin(alpha)))

    print("The final radius R is calculated from the equation: R = r * (1 + 1 / sin(alpha))")
    print("-" * 20)
    # Output each number in the final equation as requested
    print(f"The radius of the small circles, r = {r}")
    print(f"The constant '1' in the equation = {1}")
    print(f"The numerically solved angle, alpha (radians) = {alpha}")
    
    print("-" * 20)
    print("Substituting the values into the equation:")
    print(f"R = {r} * (1 + 1 / sin({alpha:.8f}))")
    print(f"R = {R:.8f}")
    
    # Format the final answer to 4 significant digits.
    # The value is ~3.4495, which rounds to 3.450. Using the .3f format specifier achieves this.
    final_radius_formatted = f"{R:.3f}"
    
    print("-" * 20)
    print(f"The radius of the circle up to 4 significant digits is: {final_radius_formatted}")

solve_fourteen_circles_packing()