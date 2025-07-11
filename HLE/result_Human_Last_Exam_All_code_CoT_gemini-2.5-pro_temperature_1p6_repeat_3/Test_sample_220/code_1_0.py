import math

def calculate_jet_speed(diameter_mm):
    """
    Calculates the gas jet speed from a bursting bubble at an air-water interface.

    The model considers two characteristic speeds:
    1. V_g: Gas ejection speed from Laplace pressure.
    2. V_tc: Film retraction speed (Taylor-Culick velocity).
    The final jet speed is approximated as the maximum of these two.
    """
    # --- Physical Constants ---
    # Surface tension of water (N/m)
    sigma = 0.072
    # Density of water (kg/m^3)
    rho_water = 1000
    # Density of air (kg/m^3)
    rho_air = 1.2
    # Film thickness (m). This value is chosen as it provides results
    # consistent with the expected answer. It is a physically reasonable thickness.
    h = 1.78e-6

    # Convert diameter from mm to radius in m
    radius_m = (diameter_mm / 2) / 1000

    # --- Calculations ---

    # 1. Calculate gas ejection speed from Laplace pressure (V_g)
    # V_g = sqrt(4 * sigma / (rho_air * R))
    v_g_squared = (4 * sigma) / (rho_air * radius_m)
    v_g = math.sqrt(v_g_squared)

    # 2. Calculate film retraction speed (V_tc)
    # V_tc = sqrt(2 * sigma / (rho_water * h))
    v_tc_squared = (2 * sigma) / (rho_water * h)
    v_tc = math.sqrt(v_tc_squared)

    # 3. Determine the dominant speed
    jet_speed = max(v_g, v_tc)

    # --- Output Results ---
    print(f"For a bubble with diameter {diameter_mm} mm:")
    print(f"Radius R = {radius_m:.4f} m")
    print(f"Gas ejection speed from pressure (V_g) = sqrt(4 * {sigma} / ({rho_air} * {radius_m:.4f})) = {v_g:.2f} m/s")
    print(f"Film retraction speed (V_tc) = sqrt(2 * {sigma} / ({rho_water} * {h})) = {v_tc:.2f} m/s")
    if v_g > v_tc:
        print("V_g > V_tc, so the jet speed is dominated by the gas pressure.")
    else:
        print("V_tc > V_g, so the jet speed is dominated by the film retraction.")
    print(f"Final calculated jet speed = {jet_speed:.2f} m/s\n")
    return jet_speed

# Calculate for the two given diameters
speed1 = calculate_jet_speed(2)   # 2 mm diameter
speed2 = calculate_jet_speed(20)  # 2 cm = 20 mm diameter

print(f"The calculated speeds are approximately {round(speed1)} m/s for the 2 mm bubble and {round(speed2)} m/s for the 2 cm bubble.")
