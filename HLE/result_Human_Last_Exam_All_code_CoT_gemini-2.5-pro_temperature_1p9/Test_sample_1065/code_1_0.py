import numpy as np
from scipy.optimize import root_scalar

def find_center_of_mass_coordinates():
    """
    This function calculates the horizontal and vertical coordinates of the center of mass
    of the string of candies.

    The physical modeling leads to a transcendental equation for the angle theta_c where
    the string lifts off the sphere: tan(theta_c) = theta_c + 0.5 * cos(theta_c).

    This code first solves this equation numerically for theta_c. Then, it uses this
    angle to calculate the coordinates of the center of mass (h, v) based on derived
    formulas which represent the weighted average of the positions of all parts of the string.
    The radius R of the pumpkin is assumed to be 1 as it cancels out from the final
    expressions for the coordinates.
    """
    
    # Define the transcendental equation to find the liftoff angle theta_c
    # f(theta) = tan(theta) - theta - 0.5 * cos(theta) = 0
    def equation_for_theta(theta):
        return np.tan(theta) - theta - 0.5 * np.cos(theta)

    # Numerically solve for theta_c in the plausible range (0, pi/2)
    # The bracket [0.1, 1.5] is a safe range for the root.
    sol = root_scalar(equation_for_theta, bracket=[0.1, 1.5])
    theta_c = sol.root

    # Calculate sin and cos for theta_c
    cos_th = np.cos(theta_c)
    sin_th = np.sin(theta_c)

    # Common denominator for the coordinate formulas
    denominator = theta_c + cos_th

    # Calculate the horizontal coordinate (h) of the center of mass
    # h = (1 - cos(theta_c) + sin(theta_c)*cos(theta_c)) / (theta_c + cos(theta_c))
    h = (1 - cos_th + sin_th * cos_th) / denominator
    
    # Calculate the vertical coordinate (v) of the center of mass
    # v = (sin(theta_c) + 0.5*cos(theta_c)**2) / (theta_c + cos(theta_c))
    v = (sin_th + 0.5 * cos_th**2) / denominator

    # The problem asks for the raw numbers of the horizontal and vertical coordinates separated by a comma.
    print(f"{h},{v}")

# Execute the function to print the result.
find_center_of_mass_coordinates()