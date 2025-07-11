import numpy as np
from scipy.integrate import quad

def analyze_poloidal_dependence():
    """
    This function demonstrates the use of a Fourier series to represent a
    function with poloidal dependence in a toroidal system.
    """

    # Let's define a sample function f(theta) that represents some
    # physical quantity's dependence on the poloidal angle theta.
    # This function is arbitrary but periodic with period 2*pi.
    def f(theta):
        # A mix of a few sinusoidal modes
        return 2.5 + 1.5 * np.cos(1 * theta) + 0.8 * np.sin(1 * theta) - 0.5 * np.cos(3 * theta)

    print("Analyzing a function's poloidal (periodic) dependence using Fourier Series.\n")
    print("The goal is to represent f(theta) as a sum of sines and cosines:")
    print("f(theta) = a0/2 + sum_{m=1 to N} [ a_m*cos(m*theta) + b_m*sin(m*theta) ]\n")

    # The number of Fourier modes to calculate
    N_modes = 5

    # Calculate the a0 coefficient (the average value)
    # a0 = (1/pi) * integral from -pi to pi of f(theta) dtheta
    # Note: The constant term in the series is a0/2
    a0_integral_result, _ = quad(f, -np.pi, np.pi)
    a0 = a0_integral_result / np.pi
    
    # Start building the equation string
    # We round for cleaner output. In real physics, full precision is used.
    equation = [f"f(theta) = {a0/2:.4f}"]

    # Calculate the am and bm coefficients for m > 0
    for m in range(1, N_modes + 1):
        # Define the functions to integrate for am and bm
        def integrand_am(theta):
            return f(theta) * np.cos(m * theta)

        def integrand_bm(theta):
            return f(theta) * np.sin(m * theta)

        # Perform the integration to find the coefficients
        # am = (1/pi) * integral of f(theta)*cos(m*theta)
        # bm = (1/pi) * integral of f(theta)*sin(m*theta)
        am_integral, _ = quad(integrand_am, -np.pi, np.pi)
        am = am_integral / np.pi

        bm_integral, _ = quad(integrand_bm, -np.pi, np.pi)
        bm = bm_integral / np.pi

        # Add non-negligible terms to our equation string
        # A small tolerance helps ignore floating point noise for zero-valued coefficients
        if abs(am) > 1e-9:
            sign = '+' if am > 0 else '-'
            equation.append(f" {sign} {abs(am):.4f}*cos({m}*theta)")
        if abs(bm) > 1e-9:
            sign = '+' if bm > 0 else '-'
            equation.append(f" {sign} {abs(bm):.4f}*sin({m}*theta)")
            
    print("Calculated Fourier Series Expansion:")
    # The ' '.join() method will print each term of the final equation
    print(' '.join(equation))

analyze_poloidal_dependence()