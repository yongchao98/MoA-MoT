import scipy.constants as const

def calculate_poynting_vector(rho, R, E_axial, v, r):
    """
    Calculates the Poynting vector for a moving, charged cylindrical rod
    in an external axial electric field.

    Args:
        rho (float): Volume charge density (C/m^3).
        R (float): Radius of the cylindrical rod (m).
        E_axial (float): Magnitude of the external axial electric field (V/m).
        v (float): Speed of the rod along its axis (m/s).
        r (float): Radial distance from the center of the rod (m).
    """
    print("-" * 50)
    print(f"Calculating Poynting vector at r = {r} m:")

    if r < 0:
        print("Radial distance r cannot be negative.")
        return

    # Epsilon_0, permeability of free space
    epsilon_0 = const.epsilon_0

    if r <= R:
        # Poynting vector inside the rod (r <= R)
        print("Point is inside or on the surface of the rod (r <= R).")
        
        # Radial component S_r = - (E * rho * v * r) / 2
        s_r_num = -(E_axial * rho * v * r)
        s_r_den = 2
        s_r = s_r_num / s_r_den
        print(f"S_r = - (E * rho * v * r) / 2")
        print(f"S_r = - ({E_axial} * {rho} * {v} * {r}) / {s_r_den}")
        print(f"S_r = {s_r:.4e} W/m^2")

        # Axial component S_z = (rho^2 * v * r^2) / (4 * epsilon_0)
        s_z_num = (rho**2 * v * r**2)
        s_z_den = (4 * epsilon_0)
        s_z = s_z_num / s_z_den
        print(f"S_z = (rho^2 * v * r^2) / (4 * epsilon_0)")
        print(f"S_z = ({rho}^2 * {v} * {r}^2) / (4 * {epsilon_0:.4e})")
        print(f"S_z = {s_z:.4e} W/m^2")

    else:
        # Poynting vector outside the rod (r > R)
        print("Point is outside the rod (r > R).")
        
        # Radial component S_r = - (E * rho * v * R^2) / (2 * r)
        s_r_num = -(E_axial * rho * v * R**2)
        s_r_den = 2 * r
        s_r = s_r_num / s_r_den
        print(f"S_r = - (E * rho * v * R^2) / (2 * r)")
        print(f"S_r = - ({E_axial} * {rho} * {v} * {R}^2) / (2 * {r})")
        print(f"S_r = {s_r:.4e} W/m^2")

        # Axial component S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
        s_z_num = (rho**2 * v * R**4)
        s_z_den = (4 * epsilon_0 * r**2)
        s_z = s_z_num / s_z_den
        print(f"S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)")
        print(f"S_z = ({rho}^2 * {v} * {R}^4) / (4 * {epsilon_0:.4e} * {r}^2)")
        print(f"S_z = {s_z:.4e} W/m^2")

    print(f"\nThe Poynting vector is S = ({s_r:.4e}) r_hat + ({s_z:.4e}) z_hat W/m^2")
    print("-" * 50 + "\n")


if __name__ == '__main__':
    # Define the physical parameters for the problem
    # Example values in SI units
    volume_charge_density = 1.0e-6  # C/m^3
    rod_radius = 0.01             # m (1 cm)
    external_E_field = 1000.0       # V/m
    rod_speed = 10.0                # m/s

    # List of radial distances to perform the calculation
    radial_distances = [0.005, 0.02] # one point inside, one outside

    # Run the calculation for each distance
    for r_val in radial_distances:
        calculate_poynting_vector(volume_charge_density, rod_radius, external_E_field, rod_speed, r_val)