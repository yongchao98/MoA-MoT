import math

def calculate_t_gate_cost():
    """
    Calculates the approximate number of physical non-Clifford (T) gates
    required to implement one fault-tolerant logical T-gate for surface
    codes of distance 3 and 5.

    It also prints the equations used for the final result.
    """
    # --- Parameters and Assumptions ---
    # Physical gate error rate p <= 1%
    P_PHYSICAL = 0.01
    # Assumed prefactor 'c' in the logical error rate formula P_L ~ c * p^((d+1)/2)
    PREFACTOR_LOGICAL_ERROR = 10.0
    # Parameters for the 15-to-1 magic state distillation protocol
    DISTILLATION_FACTOR = 35.0      # Prefactor in p_out ~ 35 * p_in^3
    DISTILLATION_INPUTS = 15        # Number of input states for the protocol

    print("This script calculates the approximate T-gate overhead for fault-tolerant computing.\n")
    print(f"Assumptions:")
    print(f"- Physical gate error rate (p): {P_PHYSICAL}")
    print(f"- Logical Clifford error rate (P_L) formula: {PREFACTOR_LOGICAL_ERROR} * p^((d+1)/2)")
    print(f"- Magic state distillation: {DISTILLATION_INPUTS}-to-1 protocol (p_out ~ {DISTILLATION_FACTOR} * p_in^3)\n")

    # --- Scenario 1: Distance-3 Code ---
    d1 = 3
    print(f"--- Scenario 1: Calculation for distance d={d1} ---")

    # 1. Calculate target logical error rate for d=3
    exponent_d1 = (d1 + 1) / 2
    target_error_d1 = PREFACTOR_LOGICAL_ERROR * (P_PHYSICAL ** exponent_d1)
    print(f"Target logical error rate for d={d1} is ~ {target_error_d1:.2e}")

    # 2. Perform distillation
    current_error_d1 = P_PHYSICAL
    cost_d1 = 1
    levels_d1 = 0
    equation_d1_parts = []
    
    # Check if the initial T-gate error is higher than the target logical error
    if current_error_d1 > target_error_d1:
        # First level of distillation
        levels_d1 += 1
        print(f"Initial T-gate error ({current_error_d1:.2e}) is too high. Applying Level {levels_d1} of distillation.")
        current_error_d1 = DISTILLATION_FACTOR * (current_error_d1 ** 3)
        cost_d1 *= DISTILLATION_INPUTS
        equation_d1_parts.append(str(DISTILLATION_INPUTS))
        print(f"New error is ~{current_error_d1:.2e}, which is below the target.")
    else:
        print("Initial T-gate error is acceptable. No distillation needed.")


    final_equation_d1 = " * ".join(equation_d1_parts)
    print("\nFor the d=3 case, the total number of non-Clifford gates is:")
    print(f"{final_equation_d1} = {cost_d1}")


    # --- Scenario 2: Distance-5 Code ---
    d2 = 5
    print(f"\n--- Scenario 2: Calculation for distance d={d2} ---")

    # 1. Calculate target logical error rate for d=5
    exponent_d2 = (d2 + 1) / 2
    target_error_d2 = PREFACTOR_LOGICAL_ERROR * (P_PHYSICAL ** exponent_d2)
    print(f"Target logical error rate for d={d2} is ~ {target_error_d2:.2e}")

    # 2. Perform distillation
    current_error_d2 = P_PHYSICAL
    cost_d2 = 1
    levels_d2 = 0
    equation_d2_parts = []
    
    # Apply levels of distillation until the error is low enough
    while current_error_d2 > target_error_d2:
        levels_d2 += 1
        print(f"Current error ({current_error_d2:.2e}) is too high. Applying Level {levels_d2} of distillation.")
        current_error_d2 = DISTILLATION_FACTOR * (current_error_d2 ** 3)
        cost_d2 *= DISTILLATION_INPUTS
        equation_d2_parts.append(str(DISTILLATION_INPUTS))
        
    final_equation_d2 = " * ".join(equation_d2_parts)
    print(f"A total of {levels_d2} levels of distillation are needed.")
    print(f"Final error is ~{current_error_d2:.2e}, which is below the target.")
    print("\nFor the d=5 case, the total number of non-Clifford gates is:")
    print(f"{final_equation_d2} = {cost_d2}")


calculate_t_gate_cost()
>>> 15 * 15 = 225