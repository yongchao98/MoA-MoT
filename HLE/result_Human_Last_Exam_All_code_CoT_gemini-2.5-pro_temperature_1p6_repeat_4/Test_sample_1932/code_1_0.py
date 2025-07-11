import math

def calculate_weight_change():
    """
    Calculates the change in weight of an hourglass based on different physical models.
    """
    # Given parameters
    d = 0.01  # diameter in meters
    h = 0.02  # height of sand column in meters
    H = 0.04  # height of each chamber in meters
    rho = 1500  # density of sand in kg/m^3
    t = 60  # time for sand to fall in seconds
    g = 9.8  # acceleration due to gravity in m/s^2

    # Calculate intermediate values
    pi = math.pi
    area = pi * d**2 / 4
    sand_mass = rho * area * h
    mass_flow_rate = sand_mass / t

    print("--- Input Parameters ---")
    print(f"Diameter (d): {d} m")
    print(f"Sand column height (h): {h} m")
    print(f"Chamber height (H): {H} m")
    print(f"Sand density (rho): {rho} kg/m^3")
    print(f"Fall time (t): {t} s")
    print("\n--- Calculated Intermediate Values ---")
    print(f"Total sand mass: {sand_mass:.6f} kg")
    print(f"Mass flow rate (ṁ): {mass_flow_rate:.6e} kg/s")

    # --- Evaluate each answer choice ---
    print("\n--- Estimated Weight Change (ΔW) for each Answer Choice ---")

    # A: ṁ * sqrt(2*g*(H-h))
    # This represents the impact force when the bottom pile has reached its full height.
    val_A_sqrt_term = math.sqrt(2 * g * (H - h))
    val_A = mass_flow_rate * val_A_sqrt_term
    print(f"A: ΔW = (π*d²*h*ρ / (4*t)) * sqrt(2*g*(H-h))")
    print(f"   Values: ({pi:.4f} * {d}² * {h} * {rho} / (4*{t})) * sqrt(2*{g}*({H}-{h}))")
    print(f"   Result: ΔW ≈ {val_A:.3e} N\n")

    # B: -ṁ * sqrt(2*g*H)
    # This represents the weight of the falling stream at the beginning of the fall.
    val_B_sqrt_term = math.sqrt(2 * g * H)
    val_B = -mass_flow_rate * val_B_sqrt_term
    print(f"B: ΔW = -(π*d²*h*ρ / (4*t)) * sqrt(2*g*H)")
    print(f"   Values: -({pi:.4f} * {d}² * {h} * {rho} / (4*{t})) * sqrt(2*{g}*{H})")
    print(f"   Result: ΔW ≈ {val_B:.3e} N\n")
    
    # C: π*d²*h²*ρ / (2*t²)
    # This represents a secondary effect related to the dynamics of the flow.
    val_C_numerator = pi * d**2 * h**2 * rho
    val_C_denominator = 2 * t**2
    val_C = val_C_numerator / val_C_denominator
    print(f"C: ΔW = (π * d² * h² * ρ) / (2 * t²)")
    print(f"   Values: ({pi:.4f} * {d}² * {h}² * {rho}) / (2 * {t}²)")
    print(f"   Result: ΔW ≈ {val_C:.3e} N\n")

    # D: ṁ * sqrt(2*g*H)
    # This represents the impact force at the very beginning of the fall.
    val_D_sqrt_term = math.sqrt(2 * g * H)
    val_D = mass_flow_rate * val_D_sqrt_term
    print(f"D: ΔW = (π*d²*h*ρ / (4*t)) * sqrt(2*g*H)")
    print(f"   Values: ({pi:.4f} * {d}² * {h} * {rho} / (4*{t})) * sqrt(2*{g}*{H})")
    print(f"   Result: ΔW ≈ {val_D:.3e} N\n")
    
    # E: 0
    # This represents the case where primary effects perfectly cancel.
    val_E = 0
    print(f"E: ΔW = 0")
    print(f"   Result: ΔW = {val_E} N\n")
    
    # F: (π*d²/4) * (H - h/2) - This isn't a force/weight change, it's a volume*length. Dimensionally incorrect.
    print(f"F: The expression is not dimensionally correct for a force.\n")

calculate_weight_change()