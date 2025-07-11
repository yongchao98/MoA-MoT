import math

def calculate_jet_speed():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface.

    The calculation considers two regimes based on the bubble size relative to the
    capillary length of water.
    """
    # --- Constants ---
    sigma = 0.072  # Surface tension of water (N/m)
    rho_air = 1.2   # Density of air (kg/m^3)
    rho_water = 1000 # Density of water (kg/m^3)
    g = 9.8         # Acceleration due to gravity (m/s^2)
    
    # --- Bubble diameters ---
    D1_mm = 2
    D2_cm = 2
    
    R1 = (D1_mm / 1000) / 2 # Radius of the first bubble in meters
    R2 = (D2_cm / 100) / 2  # Radius of the second bubble in meters

    # --- Step 1: Calculate the capillary length ---
    l_c = math.sqrt(sigma / (rho_water * g))
    print(f"Physical constants:")
    print(f"Surface tension (σ) = {sigma} N/m")
    print(f"Density of air (ρ_air) = {rho_air} kg/m^3")
    print(f"Density of water (ρ_water) = {rho_water} kg/m^3")
    print("\nFirst, we determine the capillary length (l_c) to see if bubbles are large or small.")
    print(f"l_c = sqrt(σ / (ρ_water * g)) = sqrt({sigma} / ({rho_water} * {g})) = {l_c:.4f} m\n")

    # --- Step 2: Calculate speed for the 2 mm bubble ---
    print("--- Case 1: Bubble Diameter = 2 mm ---")
    print(f"Bubble radius R1 = {R1*1000:.1f} mm. Since R1 < l_c, the bubble is small and hemispherical.")
    print(f"We use its own radius for the radius of curvature: R_c1 = {R1:.4f} m.")
    # v = sqrt(4 * σ / (ρ_air * R_c))
    v1 = math.sqrt(4 * sigma / (rho_air * R1))
    print(f"Jet speed v1 = sqrt(4 * σ / (ρ_air * R_c1)) = sqrt(4 * {sigma} / ({rho_air} * {R1})) = {v1:.2f} m/s")
    
    # --- Step 3: Calculate speed for the 2 cm bubble ---
    print("\n--- Case 2: Bubble Diameter = 2 cm ---")
    print(f"Bubble radius R2 = {R2*100:.1f} cm. Since R2 > l_c, the bubble is large and flattened by gravity.")
    print(f"We use the capillary length for the radius of curvature: R_c2 = l_c = {l_c:.4f} m.")
    # v = sqrt(4 * σ / (ρ_air * R_c))
    v2 = math.sqrt(4 * sigma / (rho_air * l_c))
    print(f"Jet speed v2 = sqrt(4 * σ / (ρ_air * R_c2)) = sqrt(4 * {sigma} / ({rho_air} * {l_c:.4f})) = {v2:.2f} m/s")

    # --- Final result ---
    print(f"\nThus, the calculated speeds are approximately {round(v1)} m/s and {round(v2)} m/s.")

calculate_jet_speed()