import scipy.constants as const

def compute_poynting_vector(R, rho, E, v, r):
    """
    Computes the Poynting vector for a moving charged cylindrical rod
    in an external electric field.

    The Poynting vector S is given in cylindrical coordinates as:
    S = S_r * r_hat + S_z * z_hat

    Args:
        R (float): Radius of the rod (m).
        rho (float): Volume charge density (C/m^3).
        E (float): Magnitude of the external axial electric field (V/m).
        v (float): Speed of the rod along its axis (m/s).
        r (float): Radial distance from the axis to compute S at (m).
    """
    print(f"--- Calculating Poynting Vector at r = {r:.3f} m ---")
    print(f"Parameters: R={R}, rho={rho}, E={E}, v={v}\n")

    if r < 0:
        print("Error: Radial distance 'r' cannot be negative.")
        return

    # The case r < R (inside the rod)
    if r < R:
        print("Calculating for a point INSIDE the rod (r < R).")
        
        # Radial component: S_r = - (E * rho * v * r) / 2
        s_r_num = - (E * rho * v * r)
        s_r_den = 2
        s_r = s_r_num / s_r_den
        
        print("Radial component (S_r):")
        print(f"  Formula: S_r = - (E * rho * v * r) / 2")
        print(f"  Calculation: S_r = - ({E} * {rho} * {v} * {r}) / {s_r_den}")
        print(f"  Result: S_r = {s_r:.4e} W/m^2\n")

        # Axial component: S_z = (rho^2 * v * r^2) / (4 * epsilon_0)
        s_z_num = (rho**2 * v * r**2)
        s_z_den = (4 * const.epsilon_0)
        s_z = s_z_num / s_z_den
        
        print("Axial component (S_z):")
        print(f"  Formula: S_z = (rho^2 * v * r^2) / (4 * epsilon_0)")
        print(f"  Calculation: S_z = ({rho:.1e}^2 * {v} * {r}^2) / (4 * {const.epsilon_0:.4e})")
        print(f"  Result: S_z = {s_z:.4e} W/m^2")

    # The case r >= R (outside the rod)
    else:
        print("Calculating for a point OUTSIDE the rod (r >= R).")
        
        # Radial component: S_r = - (E * rho * v * R^2) / (2 * r)
        s_r_num = - (E * rho * v * R**2)
        s_r_den = (2 * r)
        s_r = s_r_num / s_r_den

        print("Radial component (S_r):")
        print(f"  Formula: S_r = - (E * rho * v * R^2) / (2 * r)")
        print(f"  Calculation: S_r = - ({E} * {rho} * {v} * {R}^2) / (2 * {r})")
        print(f"  Result: S_r = {s_r:.4e} W/m^2\n")

        # Axial component: S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
        s_z_num = (rho**2 * v * R**4)
        s_z_den = (4 * const.epsilon_0 * r**2)
        s_z = s_z_num / s_z_den
        
        print("Axial component (S_z):")
        print(f"  Formula: S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)")
        print(f"  Calculation: S_z = ({rho:.1e}^2 * {v} * {R}^4) / (4 * {const.epsilon_0:.4e} * {r}^2)")
        print(f"  Result: S_z = {s_z:.4e} W/m^2")

    print("\n--- Final Result ---")
    print(f"The Poynting vector at r = {r:.3f} m is:")
    print(f"S = ({s_r:.4e}) r_hat + ({s_z:.4e}) z_hat  [W/m^2]")


# --- Example Usage ---
# You can change these values to compute the Poynting vector for different scenarios.

# Parameters for the rod and fields
rod_radius = 0.05      # R in meters
charge_density = 1e-6  # rho in Coulombs/m^3
external_E_field = 1000 # E in Volts/m
rod_speed = 100        # v in meters/s

# --- Calculation 1: Inside the rod ---
radial_distance_in = 0.025 # r in meters, where r < R
compute_poynting_vector(rod_radius, charge_density, external_E_field, rod_speed, radial_distance_in)

print("\n" + "="*60 + "\n")

# --- Calculation 2: Outside the rod ---
radial_distance_out = 0.1 # r in meters, where r > R
compute_poynting_vector(rod_radius, charge_density, external_E_field, rod_speed, radial_distance_out)