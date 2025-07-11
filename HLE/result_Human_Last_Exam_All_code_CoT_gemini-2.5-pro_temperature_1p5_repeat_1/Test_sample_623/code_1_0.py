def trace_labeled_glucose():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis and the
    Pyruvate Dehydrogenase Complex (PDC) reaction to determine the number
    of labeled CO2 molecules released.
    """
    print("--- Tracing 13C in the Catabolism of 1,4-13C Glucose ---")

    # Step 1: Define the starting glucose molecule.
    # The dictionary represents the 6 carbons of glucose, with labels at C1 and C4.
    glucose = {
        1: '13C', 2: '12C', 3: '12C',
        4: '13C', 5: '12C', 6: '12C'
    }
    print(f"\n1. Initial molecule: 1,4-13C Glucose")
    print(f"   Carbon labels: {glucose}")

    # Step 2: Glycolysis.
    # Glucose (C1-C6) -> 2 Pyruvate (3C each)
    # The mapping from glucose carbons to pyruvate carbons is specific:
    # Pyruvate's carboxyl carbon (C1) comes from Glucose's C3 or C4.
    # Pyruvate's methyl carbon (C3) comes from Glucose's C1 or C6.
    print("\n2. Glycolysis splits glucose into two 3-carbon pyruvate molecules.")
    
    # Pyruvate from top half of glucose (C1, C2, C3)
    pyruvate_1 = {'carboxyl': glucose[3], 'keto': glucose[2], 'methyl': glucose[1]}
    # Pyruvate from bottom half of glucose (C4, C5, C6)
    pyruvate_2 = {'carboxyl': glucose[4], 'keto': glucose[5], 'methyl': glucose[6]}
    
    print(f"   - Pyruvate 1 (from Glucose C1,C2,C3) has labels: {pyruvate_1}")
    print(f"   - Pyruvate 2 (from Glucose C4,C5,C6) has labels: {pyruvate_2}")

    # Step 3: Pyruvate Dehydrogenase Complex (PDC) reaction.
    # This enzyme complex removes the carboxyl group from pyruvate, releasing it as CO2.
    print("\n3. The PDC reaction converts pyruvate to acetyl-CoA, releasing the carboxyl carbon as CO2.")
    
    co2_from_pyruvate_1 = pyruvate_1['carboxyl']
    num_13co2_from_p1 = 1 if co2_from_pyruvate_1 == '13C' else 0
    print(f"   - CO2 from Pyruvate 1 (originally Glucose C3) is {co2_from_pyruvate_1}.")
    
    co2_from_pyruvate_2 = pyruvate_2['carboxyl']
    num_13co2_from_p2 = 1 if co2_from_pyruvate_2 == '13C' else 0
    print(f"   - CO2 from Pyruvate 2 (originally Glucose C4) is {co2_from_pyruvate_2}.")

    # Step 4: Citric Acid (TCA) Cycle Check
    # The remaining carbons enter the TCA cycle as acetyl-CoA.
    # The label from glucose C1 is now on the methyl group of acetyl-CoA.
    # This carbon is not released as CO2 in the first turn of the TCA cycle.
    print("\n4. The remaining labeled carbon (from Glucose C1) enters the TCA cycle but is not released as CO2 in the first turn.")

    # Step 5: Final Calculation
    print("\n--- Final Calculation ---")
    total_13co2 = num_13co2_from_p1 + num_13co2_from_p2
    print(f"Number of 13CO2 molecules from Pyruvate 1: {num_13co2_from_p1}")
    print(f"Number of 13CO2 molecules from Pyruvate 2: {num_13co2_from_p2}")
    print(f"Total 13CO2 molecules released = {num_13co2_from_p1} + {num_13co2_from_p2} = {total_13co2}")

# Execute the trace
trace_labeled_glucose()