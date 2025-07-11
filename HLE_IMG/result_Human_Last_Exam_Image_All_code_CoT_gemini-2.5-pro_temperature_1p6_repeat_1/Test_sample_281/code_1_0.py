def calculate_four_terminal_conductance():
    """
    This function prints the step-by-step derivation of the four-terminal conductance
    G_12,34 for the given quantum Hall device.
    The variables M (total edge states) and N (reflected edge states) are treated symbolically.
    """

    # Symbolic variables
    M = "M"
    N = "N"
    G0 = "e^2/h" # Conductance quantum

    print("Step 1: Model the system and edge state propagation.")
    print("-----------------------------------------------------")
    print("We use the Landauer-Büttiker formalism. The terminals are labeled 1, 5, 6, 2, 4, 3 clockwisely.")
    print("We assume clockwise edge state propagation: 1 -> 3 -> 4 -> 2 -> 6 -> 5 -> 1.")
    print(f"There are {M} spin-degenerate edge states. The QPC is located between the upper (6->5) and lower (3->4) paths.")
    print(f"At the QPC, {M-N} states are transmitted along their path, and {N} states are reflected (scattered to the opposite edge).\n")

    print("Step 2: Determine the potentials at the voltage probes.")
    print("-----------------------------------------------------")
    print("The potential of a probe is the average of the potentials of the incoming edge channels.")
    print("Source current is from 1 to 2. Let V_1 and V_2 be the potentials at contacts 1 and 2.")
    print(f"- Channels from contact 1 (at V_1) reach probe 3. So, V_3 = V_1.")
    print(f"- Channels from contact 2 (at V_2) reach probe 6. So, V_6 = V_2.")
    print(f"- Probe 4 receives {M-N} channels from probe 3 (at V_3=V_1) and {N} channels reflected from probe 6 (at V_6=V_2).")
    print(f"  V_4 = (({M}-{N})*V_3 + {N}*V_6) / {M} = (({M}-{N})*V_1 + {N}*V_2) / {M}")
    print("The measured voltage difference V_3 - V_4 is then:")
    print(f"  V_3 - V_4 = V_1 - [ (({M}-{N})*V_1 + {N}*V_2) / {M} ]")
    print(f"  V_3 - V_4 = (M*V_1 - ({M}-{N})*V_1 - {N}*V_2) / {M}")
    print(f"  V_3 - V_4 = (N*V_1 - {N}*V_2) / {M} = ({N}/{M}) * (V_1 - V_2)\n")

    print("Step 3: Determine the source current I = I_1.")
    print("-----------------------------------------------------")
    print("The current at a terminal 'k' is I_k = G0 * (M*V_k - Sum of potentials of incoming channels).")
    print("Channels arriving at terminal 1 come from probe 5. So I_1 = G0 * (M*V_1 - M*V_5).")
    print(f"- Probe 5 receives {M-N} channels from probe 6 (at V_6=V_2) and {N} channels reflected from probe 3 (at V_3=V_1).")
    print(f"  V_5 = ({N}*V_3 + ({M}-{N})*V_6) / {M} = ({N}*V_1 + ({M}-{N})*V_2) / {M}")
    print("Substituting V_5 into the equation for I_1 = I:")
    print(f"  I = ({G0}) * (M*V_1 - M*[({N}*V_1 + ({M}-{N})*V_2) / {M}])")
    print(f"  I = ({G0}) * (M*V_1 - {N}*V_1 - ({M}-{N})*V_2)")
    print(f"  I = ({G0}) * (({M}-{N})*V_1 - ({M}-{N})*V_2)")
    print(f"  I = ({G0})*({M}-{N}) * (V_1 - V_2)\n")

    print("Step 4: Calculate the four-terminal conductance G_12,34.")
    print("----------------------------------------------------------")
    print("G_12,34 = I / (V_3 - V_4)")
    print(f"G_12,34 = [ ({G0})*({M}-{N}) * (V_1 - V_2) ] / [ ({N}/{M}) * (V_1 - V_2) ]")
    print("The (V_1 - V_2) term cancels out, leaving:")
    print(f"G_12,34 = ({G0}) * ({M} * ({M} - {N})) / {N}\n")
    
    print("Final Result:")
    print("The final equation for the four-terminal conductance is:")
    final_eq = f"G_12,34 = (e^2/h) * (M * (M - N)) / N"
    print(final_eq)
    print("\nWhere the numbers and symbols in the equation represent:")
    print(f"  - The term (e^2/h) is the quantum of conductance.")
    print(f"  - The number of total edge states is {M}.")
    print(f"  - The number of transmitted states is ({M} - {N}).")
    print(f"  - The number of reflected states is {N}.")

calculate_four_terminal_conductance()