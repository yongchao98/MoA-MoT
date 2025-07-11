import numpy as np
from scipy.constants import epsilon_0

def calculate_poynting_vector(rho, E, v, R, r):
    """
    Calculates the Poynting vector for a moving cylindrical rod.

    Args:
        rho (float): Volume charge density (C/m^3).
        E (float): External electric field along the axis (V/m).
        v (float): Speed of the rod along the axis (m/s).
        R (float): Radius of the rod (m).
        r (float): Radial distance from the center of the rod (m).
    """

    print(f"Calculating Poynting vector for:\n"
          f"  Charge Density (rho) = {rho} C/m^3\n"
          f"  External E-field (E) = {E} V/m\n"
          f"  Velocity (v)         = {v} m/s\n"
          f"  Rod Radius (R)       = {R} m\n"
          f"  At radial distance (r) = {r} m\n")

    # The Poynting vector S has a radial component (S_r) and an axial component (S_z)
    # S = S_r * r_hat + S_z * k_hat

    if r < 0:
        print("Error: Radial distance 'r' cannot be negative.")
        return

    if r < R:
        # Poynting vector inside the rod
        print("Calculating for a point INSIDE the rod (r < R).")
        # S_r = - (E * rho * v * r) / 2
        s_r_coeff = - (E * rho * v * r) / 2
        
        # S_z = (rho^2 * v * r^2) / (4 * epsilon_0)
        s_z_coeff = (rho**2 * v * r**2) / (4 * epsilon_0)
        
        formula_r = "S_r = - (E * rho * v * r) / 2"
        formula_z = "S_z = (rho^2 * v * r^2) / (4 * epsilon_0)"
        
    else:
        # Poynting vector outside the rod
        print("Calculating for a point OUTSIDE the rod (r >= R).")
        # S_r = - (E * rho * v * R^2) / (2 * r)
        s_r_coeff = - (E * rho * v * R**2) / (2 * r)

        # S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
        s_z_coeff = (rho**2 * v * R**4) / (4 * epsilon_0 * r**2)
        
        formula_r = "S_r = - (E * rho * v * R^2) / (2 * r)"
        formula_z = "S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)"

    print("\n--- Formulas ---")
    print(f"Radial component: {formula_r}")
    print(f"Axial component:  {formula_z}")

    print("\n--- Final Equation ---")
    print(f"The Poynting vector at r = {r} m is:")
    # Using '+' for negative numbers to maintain the vector addition look
    print(f"S = ({s_r_coeff:.3e}) r_hat + ({s_z_coeff:.3e}) k_hat  [W/m^2]")


# --- User-defined parameters ---
# You can change these values to see how they affect the result.

# Case 1: Point inside the rod
rho_val = 1.5e-6  # C/m^3
E_val = 2000.0    # V/m
v_val = 50.0      # m/s
R_val = 0.02      # meters (2 cm)
r_val_inside = 0.01 # meters (1 cm)

calculate_poynting_vector(rho_val, E_val, v_val, R_val, r_val_inside)

print("\n" + "="*50 + "\n")

# Case 2: Point outside the rod
r_val_outside = 0.04 # meters (4 cm)
calculate_poynting_vector(rho_val, E_val, v_val, R_val, r_val_outside)