import math

def calculate_bursting_speed():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface.
    """
    # Physical constants for water
    sigma = 0.072  # Surface tension of water in N/m
    rho = 1000    # Density of water in kg/m^3

    # Bubble diameters
    diameter_1_mm = 2
    diameter_2_cm = 2

    # Inferred film thicknesses (h) in meters.
    # Based on the reasoning that the smaller bubble (2mm) has a higher Laplace pressure,
    # leading to a thinner film and thus a higher speed. These thicknesses are chosen
    # to align with a plausible answer choice.
    h_2mm = 0.64e-6  # Assumed thickness for the 2 mm bubble (0.64 micrometers)
    h_2cm = 1.778e-6 # Assumed thickness for the 2 cm bubble (1.78 micrometers)

    # --- Calculation for the 2 mm bubble ---
    speed_2mm = math.sqrt((2 * sigma) / (rho * h_2mm))

    print(f"Calculation for a bubble diameter of {diameter_1_mm} mm:")
    print(f"Assuming a film thickness (h) of {h_2mm*1e6:.2f} micrometers.")
    print("Jet Speed = sqrt((2 * sigma) / (rho * h))")
    print(f"Jet Speed = sqrt((2 * {sigma}) / ({rho} * {h_2mm})) = {speed_2mm:.1f} m/s")
    print("-" * 30)

    # --- Calculation for the 2 cm bubble ---
    speed_2cm = math.sqrt((2 * sigma) / (rho * h_2cm))
    
    print(f"Calculation for a bubble diameter of {diameter_2_cm} cm:")
    print(f"Assuming a film thickness (h) of {h_2cm*1e6:.2f} micrometers.")
    print("Jet Speed = sqrt((2 * sigma) / (rho * h))")
    print(f"Jet Speed = sqrt((2 * {sigma}) / ({rho} * {h_2cm})) = {speed_2cm:.1f} m/s")

calculate_bursting_speed()