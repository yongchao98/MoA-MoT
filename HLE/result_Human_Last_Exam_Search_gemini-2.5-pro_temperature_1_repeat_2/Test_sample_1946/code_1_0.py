import numpy as np
from scipy.integrate import quad

# --- Constants ---
g = 9.8  # m/s^2
Mc = 0.020  # Mass of cardboard core in kg
Mp = 0.200  # Mass of paper in kg
d_core = 0.04  # Diameter of core in m
r_core = d_core / 2
t = 0.0005  # Thickness of paper in m
N = 100  # Number of wraps

# --- Derived Properties ---
r_initial = r_core + N * t
A_paper = np.pi * (r_initial**2 - r_core**2)
L = A_paper / t

r_core_sq = r_core**2

# --- Helper functions that depend on unrolled length y ---
def get_radius_sq(y):
    """Calculates the square of the roll's outer radius."""
    return r_core_sq + t * (L - y) / np.pi

def get_moment_of_inertia(y, r_sq_y):
    """Calculates the moment of inertia of the roll."""
    Mp_rem = Mp * (L - y) / L
    I_core = Mc * r_core_sq # Thin shell approximation for the core
    I_paper = 0.5 * Mp_rem * (r_sq_y + r_core_sq) # Hollow cylinder
    return I_core + I_paper

# --- Integrand for the time calculation ---
def integrand(y):
    """Calculates 1/v(y), the function to be integrated."""
    if y <= 0:
        # Avoid division by zero at the start, return a large number or handle as a limit
        # The quad integrator handles this singularity
        return np.inf

    # Properties at length y
    r_sq_y = get_radius_sq(y)
    m_y = Mc + Mp * (1 - y / L)
    I_y = get_moment_of_inertia(y, r_sq_y)

    # Numerator of v^2(y) - from change in potential energy
    numerator = 2 * g * ((Mc + Mp) * y - (Mp * y**2) / (2 * L))
    
    # Denominator of v^2(y) - from kinetic energy terms
    # Includes translational and rotational KE of the roll, and KE of the unrolled sheet
    denominator = m_y + (I_y / r_sq_y) + (Mp * y) / (3 * L)
    
    if denominator <= 0: return np.inf

    v_sq = numerator / denominator
    if v_sq <= 0: return np.inf

    v = np.sqrt(v_sq)
    return 1.0 / v

# --- Perform numerical integration ---
# The quad function from scipy is suitable for this, as it can handle
# singularities at the integration limits.
time, error = quad(integrand, 0, L)

# print(f"Total length of paper (L): {L:.4f} m")
# print(f"Time to unroll: {time:.2f} seconds")
# print(f"Estimated error: {error:.2e}")

# The calculation result for time is 2.5835...