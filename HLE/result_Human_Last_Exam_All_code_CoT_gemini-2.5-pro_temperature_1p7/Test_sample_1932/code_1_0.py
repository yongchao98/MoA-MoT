import math

def solve_hourglass_weight_change():
    """
    This function calculates the estimated change in weight of a running hourglass
    based on the largest contributing effect.
    """
    # Define the parameters with given approximate values
    d = 0.01  # Diameter of the chambers in meters
    h = 0.02  # Height of the sand column when settled in meters
    H = 0.04  # Height of each chamber in meters
    rho = 1500  # Density of sand in kg/m^3
    t = 60    # Time for all the sand to fall in seconds
    g = 9.8   # Acceleration due to gravity in m/s^2

    # The largest effect influencing the weight is the initial drop as the
    # sand column is established but has not yet started to land.
    # The weight change is negative and equal to the weight of the
    # airborne sand after falling the maximum distance H.
    #
    # The formula is: Delta_W = - (mass_flow_rate) * sqrt(2 * g * H)
    # where mass_flow_rate = (total_mass) / t
    # total_mass = rho * Volume = rho * (pi * d^2 / 4) * h
    #
    # So, Delta_W = - (rho * pi * d^2 * h / (4 * t)) * sqrt(2 * g * H)
    # This corresponds to option B from the provided choices.

    # Calculate the total mass of the sand
    sand_volume = (math.pi * d**2 / 4) * h
    sand_mass = rho * sand_volume

    # Calculate the average mass flow rate
    mass_flow_rate = sand_mass / t

    # Calculate the term representing the impact velocity
    sqrt_2gH = math.sqrt(2 * g * H)

    # Calculate the change in weight
    delta_W = -mass_flow_rate * sqrt_2gH
    
    print("The formula for the estimated weight change (largest effect) is:")
    print("ΔW = - (π * d^2 * h * ρ / (4 * t)) * sqrt(2 * g * H)")
    print("\nPlugging in the given values:")
    print(f"d (diameter) = {d} m")
    print(f"h (sand height) = {h} m")
    print(f"H (chamber height) = {H} m")
    print(f"ρ (density) = {rho} kg/m^3")
    print(f"t (time) = {t} s")
    print(f"g (gravity) = {g} m/s^2")
    
    print(f"\nThe calculated weight change ΔW is: {delta_W:.4e} N")
    print("This indicates the hourglass is momentarily lighter while running.")

solve_hourglass_weight_change()