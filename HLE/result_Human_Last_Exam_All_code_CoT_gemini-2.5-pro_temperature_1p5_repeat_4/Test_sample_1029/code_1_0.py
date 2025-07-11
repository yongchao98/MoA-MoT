import numpy as np

def calculate_poynting_vector(r, R, rho, E_ext, v):
    """
    Calculates the Poynting vector S at a radial distance r from the axis.
    The function prints the calculation steps and results.

    Args:
        r (float): Radial distance from the axis (m).
        R (float): Radius of the cylindrical rod (m).
        rho (float): Volume charge density (C/m^3).
        E_ext (float): External electric field along the axis (V/m).
        v (float): Speed of the rod along the axis (m/s).
    """
    # Physical constant
    epsilon_0 = 8.854e-12  # F/m

    print("-" * 70)
    if r < R:
        print(f"Calculating Poynting vector INSIDE the rod (at r = {r} m):")
        
        # Radial component S_r
        # S_r = -E * rho * v * r / 2
        S_r = -E_ext * rho * v * r / 2
        print("S_r = -E * rho * v * r / 2")
        print(f"S_r = -({E_ext}) * ({rho}) * ({v}) * ({r}) / 2 = {S_r:.4e} W/m^2")

        # Axial component S_z
        # S_z = (rho^2 * v * r^2) / (4 * epsilon_0)
        S_z = (rho**2 * v * r**2) / (4 * epsilon_0)
        print("S_z = (rho^2 * v * r^2) / (4 * epsilon_0)")
        print(f"S_z = (({rho})^2 * ({v}) * ({r})^2) / (4 * {epsilon_0:.4e}) = {S_z:.4e} W/m^2")
        
    else: # r >= R
        print(f"Calculating Poynting vector OUTSIDE the rod (at r = {r} m):")

        # Radial component S_r
        # S_r = -E * rho * v * R^2 / (2 * r)
        S_r = -E_ext * rho * v * R**2 / (2 * r)
        print("S_r = -E * rho * v * R^2 / (2 * r)")
        print(f"S_r = -({E_ext}) * ({rho}) * ({v}) * ({R})**2 / (2 * {r}) = {S_r:.4e} W/m^2")

        # Axial component S_z
        # S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
        S_z = (rho**2 * v * R**4) / (4 * epsilon_0 * r**2)
        print("S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)")
        print(f"S_z = (({rho})^2 * ({v}) * ({R})**4) / (4 * {epsilon_0:.4e} * ({r})^2) = {S_z:.4e} W/m^2")

    print("-" * 70)


if __name__ == "__main__":
    # Define the parameters for the problem
    R_val = 0.01      # Rod radius in meters (1 cm)
    rho_val = 1e-6    # Volume charge density in C/m^3
    E_ext_val = 1000  # External electric field in V/m
    v_val = 10.0      # Speed of the rod in m/s

    print("Problem Parameters:")
    print(f"Rod Radius (R) = {R_val} m")
    print(f"Charge Density (rho) = {rho_val} C/m^3")
    print(f"External E-field (E) = {E_ext_val} V/m")
    print(f"Rod Speed (v) = {v_val} m/s")

    # --- Case 1: Point inside the rod ---
    r_inside = 0.005  # 0.5 cm
    calculate_poynting_vector(r_inside, R_val, rho_val, E_ext_val, v_val)

    # --- Case 2: Point outside the rod ---
    r_outside = 0.02  # 2 cm
    calculate_poynting_vector(r_outside, R_val, rho_val, E_ext_val, v_val)
