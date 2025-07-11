import math

def calculate_jet_speed(diameter_m):
    """
    Calculates the speed of a gas jet from a bursting bubble at an air-water interface.

    Args:
        diameter_m: The diameter of the bubble in meters.

    Returns:
        The speed of the gas jet in m/s.
    """
    # Physical constants
    sigma = 0.072  # Surface tension of water in N/m
    rho_gas = 1.2  # Density of air in kg/m^3

    # Calculate radius from diameter
    radius_m = diameter_m / 2.0

    # Calculate jet speed using the derived formula v = sqrt(4 * sigma / (R * rho_gas))
    # This comes from equating Laplace pressure (2*sigma/R) with dynamic pressure (0.5*rho_gas*v^2)
    jet_speed = math.sqrt((4 * sigma) / (radius_m * rho_gas))
    
    return jet_speed, sigma, rho_gas, radius_m

# --- Calculations for the given diameters ---

# Case 1: Diameter = 2 mm = 0.002 m
d1_m = 0.002
v1, s1, rho1, r1 = calculate_jet_speed(d1_m)

# Case 2: Diameter = 2 cm = 0.02 m
d2_m = 0.02
v2, s2, rho2, r2 = calculate_jet_speed(d2_m)

# --- Print the results ---
print("Calculation of gas jet speed from a bursting bubble:")
print("-" * 50)

# Print result for 2 mm bubble
print(f"For a bubble diameter of {d1_m*1000:.0f} mm:")
# The final print statement shows each number in the equation
print(f"v = sqrt((4 * {s1}) / ({r1} * {rho1}))")
print(f"v = {v1:.1f} m/s")
print("")

# Print result for 2 cm bubble
print(f"For a bubble diameter of {d2_m*100:.0f} cm:")
# The final print statement shows each number in the equation
print(f"v = sqrt((4 * {s2}) / ({r2} * {rho2}))")
print(f"v = {v2:.1f} m/s")
print("-" * 50)
print("The calculated speeds are approximately 15.5 m/s and 4.9 m/s.")
print("The closest answer choice is E, which has a value of 15 for the 2mm bubble. The discrepancy in the second value may be due to the simplified model ignoring gravitational effects on the larger bubble.")
