import math

def calculate_frictional_force(velocity, temperature):
    """
    Calculates a simplified frictional force in a superlubric system.
    This model demonstrates that friction increases with both velocity and temperature.

    F_friction = k * T * log(v / v0)
    where:
    - k is a proportionality constant (related to Boltzmann's constant).
    - T is the absolute temperature.
    - v is the sliding velocity.
    - v0 is a reference velocity to make the logarithm argument dimensionless.
    """
    # Define physical constants for this demonstration model
    k_proportionality = 1.0e-5  # Arbitrary units: Force/Kelvin
    v0_reference = 0.01          # Arbitrary units: m/s

    # Ensure velocity is greater than the reference velocity to avoid math domain errors
    if velocity <= v0_reference:
        return 0.0

    # Calculate the frictional force based on the model
    force = k_proportionality * temperature * math.log(velocity / v0_reference)
    return force

# --- Demonstration Case 1: Lower Temperature and Velocity ---
temp1 = 300.0  # Kelvin
vel1 = 1.0     # m/s
friction1 = calculate_frictional_force(vel1, temp1)

print("Demonstration of Frictional Response in Superlubricity:")
print("-" * 55)
# Print the equation and result for the first case
print(f"Case 1: F_friction = k * T * log(v / v0)")
print(f"Inputs: Temperature = {temp1} K, Velocity = {vel1} m/s")
print(f"Result: The calculated frictional force is {friction1:.6f} (arbitrary units).\n")

# --- Demonstration Case 2: Higher Temperature and Velocity ---
temp2 = 400.0  # Kelvin
vel2 = 10.0    # m/s
friction2 = calculate_frictional_force(vel2, temp2)

# Print the equation and result for the second case
print(f"Case 2: F_friction = k * T * log(v / v0)")
print(f"Inputs: Temperature = {temp2} K, Velocity = {vel2} m/s")
print(f"Result: The calculated frictional force is {friction2:.6f} (arbitrary units).\n")

print("Conclusion: As shown, increasing the temperature and sliding velocity leads to a higher frictional force.")
