import numpy as np

def calculate_superlubric_friction(velocity, temperature):
    """
    Calculates a simplified frictional force in a superlubric system.
    
    This model demonstrates the principle that frictional force increases with 
    both sliding velocity and temperature due to synchronized surface fluctuations.
    
    Args:
        velocity (float): The sliding velocity (e.g., in m/s).
        temperature (float): The system temperature (in Kelvin).
        
    Returns:
        float: The calculated frictional force (in arbitrary units).
    """
    # Proportionality constants for velocity and temperature
    # These are conceptual and do not represent a specific material
    k_velocity = 0.05  # Force units per m/s
    k_temperature = 0.01 # Force units per Kelvin
    
    # Calculate friction as a sum of velocity and temperature dependent terms
    friction_force = k_velocity * velocity + k_temperature * temperature
    
    # Print the equation with all the numbers
    print(f"Frictional Force = (k_v * velocity) + (k_t * temperature)")
    print(f"F = ({k_velocity} * {velocity}) + ({k_temperature} * {temperature}) = {friction_force:.4f} units")
    print("-" * 30)

def main():
    """
    Main function to run demonstrations.
    """
    print("Demonstrating how frictional force changes with velocity and temperature in a superlubric system.")
    print("Based on the principle that force increases with both factors.\n")

    # Case 1: Baseline values
    v1 = 1.0  # m/s
    t1 = 100.0 # Kelvin
    print(f"Case 1: Low Velocity ({v1} m/s), Low Temperature ({t1} K)")
    calculate_superlubric_friction(v1, t1)

    # Case 2: Increased velocity
    v2 = 5.0  # m/s
    t2 = 100.0 # Kelvin
    print(f"Case 2: High Velocity ({v2} m/s), Low Temperature ({t2} K)")
    calculate_superlubric_friction(v2, t2)

    # Case 3: Increased temperature
    v3 = 1.0  # m/s
    t3 = 300.0 # Kelvin
    print(f"Case 3: Low Velocity ({v3} m/s), High Temperature ({t3} K)")
    calculate_superlubric_friction(v3, t3)

    # Case 4: Increased both
    v4 = 5.0  # m/s
    t4 = 300.0 # Kelvin
    print(f"Case 4: High Velocity ({v4} m/s), High Temperature ({t4} K)")
    calculate_superlubric_friction(v4, t4)

if __name__ == "__main__":
    main()