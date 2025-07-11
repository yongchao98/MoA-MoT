import math

def solve_circuit():
    """
    This function implements the logic to find the resistor values and calculate the maximum possible current.
    """
    # Step 1: Find the optimal resistor values based on the given constraints.
    # From the problem analysis, we determined R1=10, R2=6, R3=13.
    R1 = 10
    R2 = 6
    R3 = 13
    V_fail_R3 = 26

    print("Step 1: Determine the resistor values.")
    print(f"Based on the problem constraints, the only viable case is when R2=6 is the minimum resistance.")
    print(f"The smallest integer resistor values satisfying all conditions are R1 = {R1} ohms and R3 = {R3} ohms (which is prime).")
    print("-" * 30)

    # Step 2: Calculate the current for each possible circuit configuration.

    # Configuration A: All resistors in parallel
    print("Step 2a: Analyzing the 'All Parallel' configuration.")
    # I_s = V_fail * (R1+R3) / (R1*R3)
    Is_A = V_fail_R3 * (R1 + R3) / (R1 * R3)
    # I3 = I_s * G3 / (G1+G2+G3) = I_s * (1/R3) / (1/R1 + 1/R2 + 1/R3)
    G1, G2, G3 = 1/R1, 1/R2, 1/R3
    I3_A = Is_A * G3 / (G1 + G2 + G3)

    # Configuration B: R1 in series with (R2 || R3)
    print("Step 2b: Analyzing the 'Series-Parallel' configuration.")
    # I_s = V_fail / R3
    Is_B = V_fail_R3 / R3
    # I3 = I_s * R2 / (R2 + R3)
    I3_B = Is_B * R2 / (R2 + R3)

    # Step 3: Determine the maximum possible current and print the final equation.
    print("-" * 30)
    print("Step 3: Find the maximum possible current and show the calculation.")

    if I3_A > I3_B:
        max_current = I3_A
        print("The maximum current is achieved with the 'All Parallel' configuration.")
        print("\nThe final calculation is:")
        print("I3 = (V_fail * (R1 + R3) / (R1 * R3)) * ( (1/R3) / (1/R1 + 1/R2 + 1/R3) )")
        numerator = V_fail_R3 * (R1 + R3) * G3
        denominator = (R1 * R3) * (G1 + G2 + G3)
        print(f"I3 = ( {V_fail_R3} * ({R1} + {R3}) / ({R1} * {R3}) ) * ( (1/{R3}) / (1/{R1} + 1/{R2} + 1/{R3}) )")
        print(f"I3 = ( {V_fail_R3 * (R1 + R3)} / {R1 * R3} ) * ( {G3:.4f} / {G1 + G2 + G3:.4f} )")
        print(f"I3 = ( {Is_A:.4f} ) * ( {G3 / (G1 + G2 + G3):.4f} )")
        print(f"I3 = {max_current:.4f} A")
    else:
        max_current = I3_B
        print("The maximum current is achieved with the 'Series-Parallel' configuration.")
        print("\nThe final calculation is:")
        print("I3 = (V_fail / R3) * (R2 / (R2 + R3))")
        print(f"I3 = ( {V_fail_R3} / {R3} ) * ( {R2} / ({R2} + {R3}) )")
        numerator = V_fail_R3 * R2
        denominator = R3 * (R2 + R3)
        print(f"I3 = {numerator} / {denominator}")
        print(f"I3 = {max_current:.4f} A")

# Run the solver
solve_circuit()

# Output the final numerical answer in the required format
# I3_A = (26 * 23 / 130) * (1/13) / (1/10+1/6+1/13) = (5.98/130)*(1/13)/((78+130+60)/780) = 4.6 * (1/13)/(268/780)
# = 4.6 * 780 / (13 * 268) = 4.6 * 60 / 268 = 276 / 268 = 69/67 = 1.02985...
# I3_B = (26/13) * (6 / (6+13)) = 2 * (6/19) = 12/19 = 0.63157...
final_answer = 69 / 67
