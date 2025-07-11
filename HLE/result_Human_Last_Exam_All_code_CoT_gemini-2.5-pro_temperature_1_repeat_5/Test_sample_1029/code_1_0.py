import scipy.constants as const

def calculate_poynting_vector(E, rho, v, R, r):
    """
    Calculates the Poynting vector components for a moving charged rod.

    Args:
        E (float): Magnitude of the external electric field (N/C).
        rho (float): Volume charge density of the rod (C/m^3).
        v (float): Speed of the rod (m/s).
        R (float): Radius of the rod (m).
        r (float): Radial distance from the axis to calculate the vector at (m).
    """
    print(f"--- Calculating Poynting Vector at r = {r} m ---")
    print(f"Parameters: E = {E} N/C, rho = {rho} C/m^3, v = {v} m/s, R = {R} m")

    if r < 0:
        print("Error: Radial distance r cannot be negative.")
        return

    if r < R:
        # Inside the rod
        print("Position is inside the rod (r < R).")
        # Radial component S_r = -E * rho * v * r / 2
        s_r = -E * rho * v * r / 2
        # Axial component S_z = rho^2 * v * r^2 / (4 * epsilon_0)
        s_z = (rho**2 * v * r**2) / (4 * const.epsilon_0)
        
        print("\nThe final equation for r < R is: S = -(E * rho * v * r / 2) r_hat + (rho^2 * v * r^2 / (4 * epsilon_0)) z_hat")
        print("Substituting the given values:")
        print(f"S_r = -({E} * {rho} * {v} * {r} / 2) = {s_r:.4e} W/m^2")
        print(f"S_z = ({rho}^2 * {v} * {r}^2 / (4 * {const.epsilon_0:.4e})) = {s_z:.4e} W/m^2")

    else:
        # Outside the rod (or on the surface)
        print("Position is on the surface or outside the rod (r >= R).")
        # Radial component S_r = -E * rho * v * R^2 / (2 * r)
        s_r = -E * rho * v * R**2 / (2 * r)
        # Axial component S_z = rho^2 * v * R^4 / (4 * epsilon_0 * r^2)
        s_z = (rho**2 * v * R**4) / (4 * const.epsilon_0 * r**2)
        
        print("\nThe final equation for r >= R is: S = -(E * rho * v * R^2 / (2 * r)) r_hat + (rho^2 * v * R^4 / (4 * epsilon_0 * r^2)) z_hat")
        print("Substituting the given values:")
        print(f"S_r = -({E} * {rho} * {v} * {R}^2 / (2 * {r})) = {s_r:.4e} W/m^2")
        print(f"S_z = ({rho}^2 * {v} * {R}^4 / (4 * {const.epsilon_0:.4e} * {r}^2)) = {s_z:.4e} W/m^2")

if __name__ == '__main__':
    # Example values for the parameters
    E_field = 100.0  # N/C
    charge_density = 0.01  # C/m^3
    speed = 10.0  # m/s
    rod_radius = 0.05  # m

    # Calculate for a point inside the rod
    radial_distance_inside = 0.025 # m
    calculate_poynting_vector(E_field, charge_density, speed, rod_radius, radial_distance_inside)
    
    print("\n" + "="*50 + "\n")

    # Calculate for a point outside the rod
    radial_distance_outside = 0.1 # m
    calculate_poynting_vector(E_field, charge_density, speed, rod_radius, radial_distance_outside)
