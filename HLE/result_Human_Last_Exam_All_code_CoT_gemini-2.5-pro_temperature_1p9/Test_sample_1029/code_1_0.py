import math

def compute_poynting_vector(rho, E_ext, v, R, r):
    """
    Computes the Poynting vector for a moving, charged cylindrical rod
    in an external electric field.

    Args:
        rho (float): Volume charge density (C/m^3)
        E_ext (float): External electric field magnitude (V/m)
        v (float): Speed of the rod (m/s)
        R (float): Radius of the rod (m)
        r (float): Radial distance from the center of the rod (m)

    Returns:
        tuple: A tuple containing the radial (S_r) and axial (S_z)
               components of the Poynting vector in W/m^2.
    """
    # Physical constants
    epsilon_0 = 8.854187817e-12  # Permittivity of free space (F/m)

    S_r = 0.0
    S_z = 0.0

    if r <= R:
        # Poynting vector inside the rod
        S_r = - (E_ext * rho * v * r) / 2
        S_z = (rho**2 * v * r**2) / (4 * epsilon_0)
    else:
        # Poynting vector outside the rod
        S_r = - (E_ext * rho * v * R**2) / (2 * r)
        S_z = (rho**2 * v * R**4) / (4 * epsilon_0 * r**2)
        
    return S_r, S_z

if __name__ == '__main__':
    # Define the parameters of the problem
    R_rod = 0.01  # Rod radius in meters (1 cm)
    rho_charge = 1.0e-6  # Volume charge density in C/m^3
    E_external = 1000.0   # External electric field in V/m
    v_rod = 10.0      # Speed of the rod in m/s

    # List of radial distances to calculate the Poynting vector at
    r_points = [0.007, 0.015] # One point inside, one point outside

    print("Calculating the Poynting vector S for the given system.\n")
    print("Parameters:")
    print(f"Rod Radius (R) = {R_rod} m")
    print(f"Volume Charge Density (ρ) = {rho_charge:.2e} C/m^3")
    print(f"External Electric Field (E) = {E_external} V/m")
    print(f"Rod Speed (v) = {v_rod} m/s\n")
    
    for r_dist in r_points:
        # Calculate the components
        s_radial, s_axial = compute_poynting_vector(rho_charge, E_external, v_rod, R_rod, r_dist)

        # Print the final equation with the computed numbers
        location = "inside" if r_dist <= R_rod else "outside"
        print(f"For r = {r_dist} m ({location} the rod):")
        # Output each number in the final equation
        print(f"Poynting Vector S = ({s_radial:.2e}) r^ + ({s_axial:.2e}) z^ [W/m^2]\n")
