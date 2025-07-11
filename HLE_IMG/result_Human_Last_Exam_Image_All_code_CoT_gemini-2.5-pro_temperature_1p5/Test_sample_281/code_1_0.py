def calculate_conductance():
    """
    This function symbolically calculates the four-terminal conductance G_12,34
    for the given quantum Hall device with a QPC.
    """
    # Print an explanation of the symbols used.
    print("This script calculates the four-terminal conductance G_12,34 symbolically.")
    print("The symbols used are:")
    print("  M: Number of spin-degenerate edge states (an integer).")
    print("  N: Number of edge states reflected by the Quantum Point Contact (QPC).")
    print("  e: The elementary charge.")
    print("  h: The Planck constant.")
    print("  G_0 = e^2/h: The quantum of conductance.\n")

    # The derivation steps are based on the Landauer-Büttiker formalism.
    print("--- Derivation Steps ---")

    # Step 1: Define total channels
    print("1. The total number of edge channels, M_total, is 2*M due to spin degeneracy.")
    print("   M_total = 2 * M\n")

    # Step 2: Calculate V_34
    print("2. The potential at the voltage probes (3 and 4) is determined by the potentials of the incoming edge channels.")
    print("   Let V_1 and V_2 be the potentials at the current source (terminal 1) and drain (terminal 2).")
    print("   Edge channels from terminal 2 arrive at probe 4, so V_4 = V_2.")
    print("   Probe 3 receives N reflected channels from the V_1 side and (M_total - N) transmitted channels from the V_2 side.")
    print("   Thus, V_3 = (N * V_1 + (2*M - N) * V_2) / (2*M)")
    print("   The voltage difference V_34 = V_3 - V_4 is therefore:")
    print("   V_34 = (N / (2 * M)) * (V_1 - V_2)\n")

    # Step 3: Calculate Current I
    print("3. The current I flowing from terminal 1 to 2 is determined by the number of channels transmitted through the device.")
    print("   The net current is carried by the (M_total - N) channels that are NOT reflected by the QPC.")
    print("   I = (M_total - N) * G_0 * (V_1 - V_2)")
    print("   I = (2*M - N) * G_0 * (V_1 - V_2)\n")

    # Step 4: Calculate the final conductance G_12,34
    print("4. The four-terminal conductance G_12,34 is the ratio I / V_34.")
    print("   G_12,34 = [ (2*M - N) * G_0 * (V_1 - V_2) ] / [ (N / (2 * M)) * (V_1 - V_2) ]")
    print("   The (V_1 - V_2) terms cancel out.\n")

    # Final Result
    print("--- Final Formula ---")
    final_formula = "G_12,34 = (e^2/h) * (2*M * (2*M - N)) / N"
    print(final_formula)
    print("")

    # Breakdown as requested
    print("--- Breakdown of the Final Equation ---")
    print("The final equation has the form G = (e^2/h) * C, where C is a coefficient given by a fraction:")
    print("  Numerator = term_A * (term_B - term_C)")
    print("  Denominator = term_C\n")
    print("The terms in the equation are composed of the given variables M and N, and the number 2:")
    print("  term_A represents the total number of channels: 2 * M")
    print("  term_B represents the total number of channels: 2 * M")
    print("  term_C represents the number of reflected channels: N")

calculate_conductance()