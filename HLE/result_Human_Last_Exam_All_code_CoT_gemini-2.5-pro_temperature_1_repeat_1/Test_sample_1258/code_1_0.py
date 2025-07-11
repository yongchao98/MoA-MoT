import numpy as np
from scipy.special import ellipk, ellipe

def calculate_demagnetizing_factor(g):
    """
    Calculates and displays the fluxmetric demagnetizing factor for a magnetic cylinder.

    The analytical expression for the fluxmetric demagnetizing factor (N_f) for a
    uniformly magnetized cylinder is:
    N_f = (4 / (g * pi)) * [K(k) - E(k)]
    where g is the length-to-diameter ratio, K and E are complete elliptic
    integrals of the first and second kind, and the modulus k is given by:
    k^2 = 1 / (1 + g^2 / 4)

    Args:
        g (float): The length-to-diameter ratio (L/D).
    """
    if g <= 0:
        print("Error: Length-to-diameter ratio 'g' must be a positive number.")
        return

    # Step 1: Calculate k^2 from g
    # This is the parameter 'm' for scipy's elliptic integral functions.
    k_squared = 1.0 / (1.0 + g**2 / 4.0)

    # Step 2: Calculate the complete elliptic integrals K(k) and E(k)
    # scipy.special.ellipk(m) computes K(k) where m = k^2
    # scipy.special.ellipe(m) computes E(k) where m = k^2
    K_val = ellipk(k_squared)
    E_val = ellipe(k_squared)

    # Step 3: Calculate the fluxmetric demagnetizing factor N_f
    N_f = (4.0 / (g * np.pi)) * (K_val - E_val)

    # Print the detailed breakdown of the calculation
    print("The analytical expression for the fluxmetric demagnetizing factor is:")
    print("N_f = (4 / (g * pi)) * [K(k) - E(k)], with k^2 = 1 / (1 + g^2/4)")
    print("\n--- Calculation for g = {} ---".format(g))
    
    # Print the calculation of k^2
    print("1. Calculate k^2:")
    print("k^2 = 1 / (1 + {:.4f}^2 / 4) = {:.4f}".format(g, k_squared))

    # Print the values of the elliptic integrals
    print("\n2. Calculate elliptic integrals K(k) and E(k) for k^2 = {:.4f}:".format(k_squared))
    print("K(k) = {:.4f}".format(K_val))
    print("E(k) = {:.4f}".format(E_val))
    
    # Print the final equation with all the numbers
    print("\n3. Calculate N_f:")
    print("N_f = (4 / ({:.4f} * {:.4f})) * ({:.4f} - {:.4f})".format(g, np.pi, K_val, E_val))
    
    # Print the final result
    print("\nFinal Result:")
    print("N_f = {:.4f}".format(N_f))


# --- User Input ---
# You can change the value of g_ratio to calculate the factor for a different geometry.
g_ratio = 2.0

# Execute the calculation
calculate_demagnetizing_factor(g_ratio)