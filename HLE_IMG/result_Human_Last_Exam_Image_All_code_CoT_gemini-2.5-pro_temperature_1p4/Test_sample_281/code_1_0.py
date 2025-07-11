def solve_conductance():
    """
    Calculates and explains the derivation of the four-terminal conductance
    G_12,34 for the given quantum Hall device.
    """
    # Define symbolic variables as strings for printing
    M = 'M'
    N = 'N'
    V1 = 'V1'
    V2 = 'V2'
    G0 = 'e^2/h'

    print("This script calculates the four-terminal conductance G_12,34 of the given quantum Hall device.")
    print("The calculation is based on the Landauer-Büttiker formalism for edge state transport.")
    print(f"Let M be the total number of spin-degenerate edge states.")
    print(f"Let N be the number of edge states reflected by the Quantum Point Contact (QPC).")
    print(f"The number of transmitted states is therefore ({M} - {N}).")
    print(f"The fundamental unit of conductance is the conductance quantum, G_0 = {G0}.")
    print("-" * 50)

    print("Step 1: Determine the voltages at the probes.")
    print("Current is sourced at terminal 1 (voltage V1) and drained at terminal 2 (voltage V2).")
    print("Edge states travel clockwise: 1 -> 5 -> QPC -> 6 -> 2 -> 4 -> QPC -> 3 -> 1.")
    print("Voltage probes (terminals 3, 4, 5, 6) draw no current, so their voltage equilibrates with incoming channels.")
    print(f"V5 = {V1} (receives {M} channels from terminal 1)")
    print(f"V6 = {V1} (receives {M}-{N} transmitted channels from the stream of terminal 1)")
    print(f"V4 = {V2} (receives {M}-{N} channels from terminal 2)")
    print("-" * 50)

    print("Step 2: Calculate V3.")
    print("Terminal 3 receives two sets of channels:")
    print(f"- {N} channels reflected by the QPC from the {V1} stream.")
    print(f"- {M}-{N} channels coming from terminal 4 (at {V2}).")
    print("The voltage V3 is the weighted average of the incoming potentials:")
    print(f"V3 = ({N}*{V1} + ({M}-{N})*{V2}) / ({N} + ({M}-{N}))")
    print(f"V3 = ({N}*{V1} + ({M}-{N})*{V2}) / {M}")
    print("-" * 50)
    
    print("Step 3: Calculate the measured voltage difference V_34 = V3 - V4.")
    print(f"V_34 = V3 - {V2}")
    print(f"V_34 = [({N}*{V1} + ({M}-{N})*{V2}) / {M}] - {V2}")
    print(f"V_34 = ({N}*{V1} + ({M}-{N})*{V2} - {M}*{V2}) / {M}")
    print(f"V_34 = ({N}*{V1} - {N}*{V2}) / {M}")
    print(f"V_34 = ({N}/{M}) * ({V1} - {V2})")
    print("-" * 50)

    print("Step 4: Calculate the total current I.")
    print("The current I flows from terminal 1 to 2. We can calculate it based on the net flow into terminal 2.")
    print(f"I = (number of channels from 1 to 2) * G_0 * ({V1} - {V2})")
    print(f"The number of channels connecting the path from 1 to 2 is the number transmitted through the QPC, {M} - {N}.")
    print(f"I = ({M}-{N}) * ({G0}) * ({V1} - {V2})")
    print("-" * 50)
    
    print("Step 5: Calculate the conductance G_12,34 = I / V_34.")
    print(f"G_12,34 = [({M}-{N}) * ({G0}) * ({V1} - {V2})] / [({N}/{M}) * ({V1} - {V2})]")
    print(f"The ({V1} - {V2}) terms cancel out.")
    print(f"G_12,34 = ({G0}) * ({M}-{N}) / ({N}/{M})")
    print("-" * 50)
    
    print("Final result:")
    print("The four-terminal conductance is G_12,34 = (e^2/h) * [M * (M - N)] / N")
    print("\nDecomposition of the final equation's terms:")
    numerator_part_1 = M
    numerator_part_2 = M
    numerator_part_3 = N
    denominator_part_1 = N
    print(f"G_12,34 = (e^2/h) * (Term1 * (Term2 - Term3)) / (Term4)")
    print(f"Term1 = {numerator_part_1}")
    print(f"Term2 = {numerator_part_2}")
    print(f"Term3 = {numerator_part_3}")
    print(f"Term4 = {denominator_part_1}")
    print("\nFinal equation with each term shown:")
    print(f"G_12,34 = ({G0}) * ( {numerator_part_1} * ( {numerator_part_2} - {numerator_part_3} ) ) / {denominator_part_1}")

solve_conductance()