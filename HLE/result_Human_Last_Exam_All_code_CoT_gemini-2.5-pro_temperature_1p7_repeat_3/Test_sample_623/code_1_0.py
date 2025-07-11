def trace_glycolysis_labels():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and pyruvate decarboxylation
    to determine the number of labeled CO2 molecules released.
    """
    print("Step 1: Define the starting molecule, 1,4-13C glucose.")
    # In 1,4-13C glucose, carbons at position 1 and 4 are labeled.
    # We represent a labeled carbon with '[13C]' and an unlabeled one with '[C]'.
    # Glucose carbons: C1-C2-C3-C4-C5-C6
    glucose = {1: '[13C]', 2: '[C]', 3: '[C]', 4: '[13C]', 5: '[C]', 6: '[C]'}
    print(f"Starting Glucose: {glucose[1]}-{glucose[2]}-{glucose[3]}-{glucose[4]}-{glucose[5]}-{glucose[6]}\n")

    print("Step 2: Glycolysis - Cleavage of Fructose-1,6-bisphosphate.")
    # Aldolase cleaves the 6-carbon sugar into two 3-carbon molecules:
    # - Dihydroxyacetone phosphate (DHAP) from carbons 1, 2, 3 of glucose.
    # - Glyceraldehyde-3-phosphate (GAP) from carbons 4, 5, 6 of glucose.
    dhap = {1: glucose[1], 2: glucose[2], 3: glucose[3]}
    gap_1 = {1: glucose[4], 2: glucose[5], 3: glucose[6]}
    print(f"DHAP (from C1,C2,C3):   {dhap[1]}-{dhap[2]}-{dhap[3]}")
    print(f"GAP (from C4,C5,C6):     {gap_1[1]}-{gap_1[2]}-{gap_1[3]}\n")

    print("Step 3: Glycolysis - Isomerization of DHAP to GAP.")
    # Triose phosphate isomerase converts DHAP into a second molecule of GAP.
    # The carbon positions map directly: DHAP C1 -> GAP C1, etc.
    gap_2 = dhap
    print(f"We now have two molecules of Glyceraldehyde-3-Phosphate (GAP):")
    print(f" - GAP Molecule 1: {gap_1[1]}-{gap_1[2]}-{gap_1[3]}")
    print(f" - GAP Molecule 2 (from DHAP): {gap_2[1]}-{gap_2[2]}-{gap_2[3]}\n")

    print("Step 4: Glycolysis - Conversion of GAP to Pyruvate.")
    # Each GAP molecule is converted to a pyruvate molecule.
    # The C1 of GAP (aldehyde) becomes the C1 of pyruvate (carboxyl group).
    pyruvate_1 = gap_1
    pyruvate_2 = gap_2
    print(f"This results in two molecules of Pyruvate:")
    print(f" - Pyruvate 1 (from GAP 1): Carbons C1,C2,C3 are {pyruvate_1[1]}, {pyruvate_1[2]}, {pyruvate_1[3]}")
    print(f" - Pyruvate 2 (from GAP 2): Carbons C1,C2,C3 are {pyruvate_2[1]}, {pyruvate_2[2]}, {pyruvate_2[3]}\n")

    print("Step 5: Pyruvate Decarboxylation (Link Reaction).")
    # This is the first step where CO2 is released.
    # The C1 carbon (carboxyl group) of pyruvate is removed and released as CO2.
    labeled_co2_count = 0
    
    print("Analyzing Pyruvate Molecule 1...")
    co2_from_pyruvate_1 = pyruvate_1[1]
    if co2_from_pyruvate_1 == '[13C]':
        labeled_co2_count += 1
    print(f"Pyruvate 1 C1 is {co2_from_pyruvate_1}. CO2 released is {('labeled' if co2_from_pyruvate_1 == '[13C]' else 'unlabeled')}.")

    print("Analyzing Pyruvate Molecule 2...")
    co2_from_pyruvate_2 = pyruvate_2[1]
    if co2_from_pyruvate_2 == '[13C]':
        labeled_co2_count += 1
    print(f"Pyruvate 2 C1 is {co2_from_pyruvate_2}. CO2 released is {('labeled' if co2_from_pyruvate_2 == '[13C]' else 'unlabeled')}.\n")
    
    unlabeled_co2_count = 2 - labeled_co2_count

    print("--- Final Calculation ---")
    print(f"Each of the 2 pyruvate molecules produced from one glucose molecule has its C1 carbon labeled with 13C.")
    print(f"The C1 carbon is lost as CO2 during pyruvate decarboxylation.")
    print(f"Therefore, the total number of 13C-labeled CO2 molecules released is {labeled_co2_count}.")
    print("\nFinal Equation:")
    print(f"1 [1,4-13C Glucose] -> 2 [1-13C Pyruvate] -> {labeled_co2_count} [13CO2] + {unlabeled_co2_count} [CO2]")


trace_glycolysis_labels()
<<<2>>>