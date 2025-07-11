import numpy as np

def calculate_fourier_series_demo():
    """
    This function demonstrates the Fourier series expansion for a sample
    periodic function representing poloidal dependence.
    
    A function f(theta) on the interval [0, 2*pi] can be represented as:
    f(theta) = a0/2 + sum_{m=1 to inf} [am * cos(m*theta) + bm * sin(m*theta)]
    
    We will analyze the function: f(theta) = 2.0 * sin(1*theta) + 0.5 * cos(3*theta)
    The non-zero coefficients should be b1 = 2.0 and a3 = 0.5.
    
    The coefficients are calculated by the integrals:
    am = (1/pi) * integral from 0 to 2*pi of f(theta)*cos(m*theta) d(theta)
    bm = (1/pi) * integral from 0 to 2*pi of f(theta)*sin(m*theta) d(theta)
    """
    
    # Define the periodic function f(theta)
    def f(theta):
        return 2.0 * np.sin(1 * theta) + 0.5 * np.cos(3 * theta)

    # Use numerical integration (trapezoidal rule) to find coefficients
    N_points = 1000  # Number of points for integration
    theta = np.linspace(0, 2 * np.pi, N_points, endpoint=False)
    d_theta = theta[1] - theta[0]

    # Calculate a0
    a0 = (1 / np.pi) * np.sum(f(theta) * d_theta)

    # We expect an equation of the form:
    # f(theta) ≈ a0/2 + a1*cos(theta) + b1*sin(theta) + a2*cos(2*theta) + b2*sin(2*theta) + ...
    # Let's find the first few coefficients.
    
    print("Calculating Fourier coefficients (am, bm) for f(theta) = 2.0*sin(theta) + 0.5*cos(3*theta)")
    print("-" * 70)
    print("Final Equation Form: f(theta) = a0/2 + sum[ am*cos(m*theta) + bm*sin(m*theta) ]\n")
    
    print(f"a0 = {a0:.5f}")
    
    max_m = 5 # Calculate up to the 5th mode
    a_coeffs = []
    b_coeffs = []
    
    for m in range(1, max_m + 1):
        # Calculate am
        integrand_a = f(theta) * np.cos(m * theta)
        am = (1 / np.pi) * np.sum(integrand_a * d_theta)
        
        # Calculate bm
        integrand_b = f(theta) * np.sin(m * theta)
        bm = (1 / np.pi) * np.sum(integrand_b * d_theta)
        
        a_coeffs.append(am)
        b_coeffs.append(bm)
        
        # Output each number in the final equation's components
        print(f"For mode m={m}:")
        print(f"  Component: {am:+.5f} * cos({m}*theta) + {bm:+.5f} * sin({m}*theta)")

calculate_fourier_series_demo()