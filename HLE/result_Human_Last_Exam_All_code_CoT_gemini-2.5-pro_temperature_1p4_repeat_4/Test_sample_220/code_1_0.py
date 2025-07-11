import math

def calculate_jet_speed():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface.

    The calculation is based on the Taylor-Culick velocity for film retraction,
    assuming the film thickness 'h' scales with the square root of the bubble diameter 'D'.
    """

    # Define physical constants for the air-water system
    sigma = 0.072  # Surface tension of water (N/m)
    rho = 1000     # Density of water (kg/m^3)
    
    # Proportionality constant for film thickness model: h = k * sqrt(D)
    # This empirical value is chosen based on experimental data for such systems.
    k = 1.34e-5    # Units: m^0.5

    # Given bubble diameters
    D1_mm = 2
    D2_cm = 2

    # Convert diameters to SI units (meters)
    D1_m = D1_mm / 1000.0
    D2_m = D2_cm / 100.0

    print("Calculating jet speed for a bursting bubble with diameter of 2 mm:")
    
    # --- Calculation for D = 2 mm ---
    # 1. Calculate film thickness h1 using the scaling law
    h1 = k * math.sqrt(D1_m)
    
    # 2. Calculate jet speed V1 using the Taylor-Culick formula
    V1 = math.sqrt((2 * sigma) / (rho * h1))

    print(f"Film thickness (h) = {k:.2e} * sqrt({D1_m}) = {h1:.3e} m")
    print(f"Jet Speed (V) = sqrt(2 * {sigma} / ({rho} * {h1:.3e})) = {V1:.1f} m/s")
    print("-" * 40)

    print("Calculating jet speed for a bursting bubble with diameter of 2 cm:")

    # --- Calculation for D = 2 cm ---
    # 1. Calculate film thickness h2 using the scaling law
    h2 = k * math.sqrt(D2_m)

    # 2. Calculate jet speed V2 using the Taylor-Culick formula
    V2 = math.sqrt((2 * sigma) / (rho * h2))
    
    print(f"Film thickness (h) = {k:.2e} * sqrt({D2_m}) = {h2:.3e} m")
    print(f"Jet Speed (V) = sqrt(2 * {sigma} / ({rho} * {h2:.3e})) = {V2:.1f} m/s")

calculate_jet_speed()