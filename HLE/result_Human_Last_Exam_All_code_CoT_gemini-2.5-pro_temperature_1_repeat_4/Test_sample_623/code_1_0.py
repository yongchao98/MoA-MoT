def trace_glucose_metabolism():
    """
    Traces the labeled carbons from 1,4-13C glucose through glycolysis and
    pyruvate decarboxylation to determine the number of labeled CO2 molecules released.
    """
    # Step 1: Define the labeled 1,4-13C glucose molecule.
    # '13C' indicates a labeled carbon, '12C' is the standard isotope.
    glucose_carbons = {
        1: '13C', 2: '12C', 3: '12C',
        4: '13C', 5: '12C', 6: '12C'
    }
    print("Step 1: The initial molecule is 1,4-13C glucose.")
    print(f"Glucose carbon labels: {glucose_carbons}\n")

    # Step 2: Simulate the cleavage of glucose in glycolysis.
    # Glucose is split into two 3-carbon molecules (triose phosphates).
    # One gets C1, C2, C3 and the other gets C4, C5, C6.
    # These are both converted to pyruvate.
    print("Step 2: Glycolysis splits glucose into two pyruvate molecules.")
    
    # Pyruvate from the top half (C1, C2, C3) of glucose
    # Glucose C1 becomes the carboxyl carbon of the first pyruvate.
    pyruvate1_carboxyl_source = 1
    pyruvate1_carboxyl_label = glucose_carbons[pyruvate1_carboxyl_source]
    print(f"The first pyruvate's carboxyl group comes from glucose C{pyruvate1_carboxyl_source}, which is labeled '{pyruvate1_carboxyl_label}'.")

    # Pyruvate from the bottom half (C4, C5, C6) of glucose
    # Glucose C4 becomes the carboxyl carbon of the second pyruvate.
    pyruvate2_carboxyl_source = 4
    pyruvate2_carboxyl_label = glucose_carbons[pyruvate2_carboxyl_source]
    print(f"The second pyruvate's carboxyl group comes from glucose C{pyruvate2_carboxyl_source}, which is labeled '{pyruvate2_carboxyl_label}'.\n")
    
    # Step 3: Simulate CO2 release from pyruvate decarboxylation.
    # The carboxyl group of each pyruvate is released as CO2.
    print("Step 3: Pyruvate decarboxylation releases the carboxyl group as CO2.")
    co2_from_pyruvate1 = pyruvate1_carboxyl_label
    co2_from_pyruvate2 = pyruvate2_carboxyl_label
    print(f"CO2 from the first pyruvate is therefore '{co2_from_pyruvate1}'.")
    print(f"CO2 from the second pyruvate is therefore '{co2_from_pyruvate2}'.\n")

    # Step 4: Count the number of 13C-labeled CO2 molecules.
    print("Step 4: Count the total number of 13CO2 molecules.")
    labeled_co2_count = 0
    contribution_from_c1 = 0
    contribution_from_c4 = 0

    if co2_from_pyruvate1 == '13C':
        labeled_co2_count += 1
        contribution_from_c1 = 1
    
    if co2_from_pyruvate2 == '13C':
        labeled_co2_count += 1
        contribution_from_c4 = 1

    # Final equation and result
    print("Final Equation:")
    print(f"Labeled CO2 from glucose C1 origin: {contribution_from_c1}")
    print(f"Labeled CO2 from glucose C4 origin: {contribution_from_c4}")
    print(f"Total labeled 13CO2 molecules = {contribution_from_c1} + {contribution_from_c4} = {labeled_co2_count}")

# Execute the function to get the answer.
trace_glucose_metabolism()