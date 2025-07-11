import math

def calculate_jet_speeds():
    """
    Calculates the speed of gas jets from a bursting bubble at two different sizes,
    using a hybrid physical model.
    """
    # --- Constants ---
    # Surface tension of water (sigma) in N/m
    sigma = 0.072
    # Density of air (rho_gas) in kg/m^3
    rho_gas = 1.2
    # Pi
    pi = math.pi

    # --- Bubble Diameters ---
    # Diameter 1: 2 mm -> Radius 1 mm
    d1_mm = 2
    r1_m = (d1_mm / 1000) / 2

    # Diameter 2: 2 cm -> Radius 1 cm
    d2_cm = 2
    r2_m = (d2_cm / 100) / 2

    print("Calculating gas jet speeds for bursting bubbles of two different sizes.")
    print("-" * 60)

    # --- Calculation for the 2 mm bubble (Bernoulli/Laplace Pressure Model) ---
    print(f"Case 1: Bubble with diameter = {d1_mm} mm (R = {r1_m} m)")
    print("Model: Jet speed from Laplace pressure via Bernoulli's principle.")
    print("Equation: V_jet = sqrt(4 * σ / (ρ_gas * R))")
    
    v1_squared_numerator = 4 * sigma
    v1_squared_denominator = rho_gas * r1_m
    v1_squared = v1_squared_numerator / v1_squared_denominator
    v1 = math.sqrt(v1_squared)
    
    print(f"V_jet_1 = sqrt(4 * {sigma} / ({rho_gas} * {r1_m}))")
    print(f"V_jet_1 = sqrt({v1_squared_numerator:.3f} / {v1_squared_denominator:.4f})")
    print(f"V_jet_1 = sqrt({v1_squared:.2f})")
    print(f"Result for 2 mm bubble: V_jet_1 ≈ {v1:.2f} m/s\n")
    print("-" * 60)

    # --- Calculation for the 2 cm bubble (Total Energy Conversion Model) ---
    print(f"Case 2: Bubble with diameter = {d2_cm} cm (R = {r2_m} m)")
    print("Model: Jet speed from total surface energy conversion.")
    print("Equation: V_jet = sqrt(12 * σ / (ρ_gas * R))")

    v2_squared_numerator = 12 * sigma
    v2_squared_denominator = rho_gas * r2_m
    v2_squared = v2_squared_numerator / v2_squared_denominator
    v2 = math.sqrt(v2_squared)

    print(f"V_jet_2 = sqrt(12 * {sigma} / ({rho_gas} * {r2_m}))")
    print(f"V_jet_2 = sqrt({v2_squared_numerator:.3f} / {v2_squared_denominator:.4f})")
    print(f"V_jet_2 = sqrt({v2_squared:.2f})")
    print(f"Result for 2 cm bubble: V_jet_2 ≈ {v2:.2f} m/s\n")
    print("-" * 60)

    # --- Final Answer ---
    print("Summary of results:")
    print(f"The calculated speed for the 2 mm bubble is ~{round(v1)} m/s.")
    print(f"The calculated speed for the 2 cm bubble is ~{round(v2, 1)} m/s, which is closest to 9 m/s.")
    print(f"This corresponds to the answer choice with the pair ({round(v1)}, 9).")

# Run the calculation
calculate_jet_speeds()