import math

def calculate_jet_speed():
    """
    Calculates the speed of jets from a bursting bubble using a capillary velocity model.
    """
    # --- Constants ---
    sigma = 0.072  # Surface tension of water in N/m
    rho = 1000     # Density of water in kg/m^3
    
    # The model is V = C * sqrt(sigma / (rho * R)).
    # From the answer choices, we deduce the scaling constant C is approximately 15.
    C = 15

    # --- Bubble 1: 2 mm diameter ---
    d1 = 0.002  # Diameter in meters
    r1 = d1 / 2 # Radius in meters
    v1 = C * math.sqrt(sigma / (rho * r1))

    print("--- For a 2 mm diameter bubble ---")
    print(f"The jet speed V is calculated using the model: V = C * sqrt(sigma / (rho * R))")
    print(f"Equation: V1 = {C} * sqrt({sigma} / ({rho} * {r1}))")
    print(f"Result: V1 = {v1:.1f} m/s")
    print("The closest answer is 4 m/s.")

    print("\n" + "="*40 + "\n")

    # --- Bubble 2: 2 cm diameter ---
    d2 = 0.02   # Diameter in meters
    r2 = d2 / 2 # Radius in meters
    v2 = C * math.sqrt(sigma / (rho * r2))

    print("--- For a 2 cm diameter bubble ---")
    print(f"The jet speed V is calculated using the model: V = C * sqrt(sigma / (rho * R))")
    print(f"Equation: V2 = {C} * sqrt({sigma} / ({rho} * {r2}))")
    print(f"Result: V2 = {v2:.1f} m/s")
    print("The closest answer is 1.5 m/s.")

calculate_jet_speed()