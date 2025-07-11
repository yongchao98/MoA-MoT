import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit

def analyze_integral_growth():
    """
    This function numerically analyzes the growth of the integral of |nabla u|^2
    for the 1D solution u(x1) = tanh(x1/sqrt(2)) to the Allen-Cahn equation.
    It computes the integral for various radii R and fits the result to a
    power law to find the growth exponent.
    """

    # The squared magnitude of the gradient for u(x1) = tanh(x1/sqrt(2))
    # |nabla u|^2 = 0.5 * sech(x1/sqrt(2))^4
    # The integral over B_R is Integrate[pi*(R^2-x1^2)*|nabla u|^2, {x1, -R, R}]
    def integrand(x1, R):
        sech_val = 1 / np.cosh(x1 / np.sqrt(2))
        return np.pi * (R**2 - x1**2) * 0.5 * sech_val**4

    # Radii to test
    R_values = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    integral_values = []

    print("Calculating integral for different radii R:")
    for R in R_values:
        val, err = quad(integrand, -R, R, args=(R,))
        integral_values.append(val)
        print(f"R = {R:.1f}, Integral = {val:.4f}")

    integral_values = np.array(integral_values)

    # Fit the data to a power law: I(R) = C * R^a
    # This is equivalent to fitting a line to log(I) vs log(R)
    # log(I) = log(C) + a * log(R)
    
    # Define the power law function for fitting
    def power_law(x, a, c):
        return c * np.power(x, a)
        
    popt, pcov = curve_fit(power_law, R_values, integral_values, p0=[2.0, 1.0])
    
    exponent_a = popt[0]

    print("\nFitting integral values to the curve C * R^a ...")
    print(f"The calculated exponent 'a' for the 1D solution is: {exponent_a:.4f}")
    
    print("\nThis numerical result supports the theoretical analysis that for the 1D solution,")
    print("the integral grows like R^2.")
    print("However, as explained in the text, there exist other solutions (e.g., radial solutions)")
    print("for which the integral grows only like R^1. The existence of these solutions")
    print("constrains the value of 'a' to be at most 1.")
    print("\nTherefore, the largest possible value for 'a' is 1.")

analyze_integral_growth()
