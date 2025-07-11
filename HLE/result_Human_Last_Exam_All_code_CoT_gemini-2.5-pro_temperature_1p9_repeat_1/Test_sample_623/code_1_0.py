def trace_glycolysis_co2_release():
    """
    This script traces the 13C labels from 1,4-13C glucose through glycolysis
    and subsequent decarboxylation to determine the number of labeled CO2 molecules released.
    """
    # Step 1: Define the 1,4-13C glucose molecule.
    # We represent the 6 carbons of glucose in a list. C1 is at index 0.
    # 13 indicates a 13C label, 12 indicates the standard 12C.
    glucose = [12] * 6
    glucose[0] = 13  # Label at C1
    glucose[3] = 13  # Label at C4

    print("--- Step 1: Starting Molecule ---")
    glucose_str = "-".join(["13C" if c == 13 else "C" for c in glucose])
    print(f"The carbon chain of 1,4-13C glucose is represented as: C1-C2-C3-C4-C5-C6")
    print(f"Labeled Glucose: {glucose_str}\n")

    # Step 2: Trace carbons through glycolysis.
    # Glycolysis cleaves glucose (C1-C2-C3-C4-C5-C6) between C3 and C4.
    # This creates two 3-carbon molecules that are converted to pyruvate.
    # - Pyruvate A originates from glucose carbons C1, C2, and C3.
    # - Pyruvate B originates from glucose carbons C4, C5, and C6.
    # The carboxyl carbon of Pyruvate (which becomes CO2) comes from:
    # - C1 of glucose for Pyruvate A.
    # - C4 of glucose for Pyruvate B.

    pyruvate_A_carboxyl_label = glucose[0]
    pyruvate_B_carboxyl_label = glucose[3]
    
    print("--- Step 2: Glycolysis ---")
    print("Glucose is split into two 3-carbon pyruvate molecules.")
    print("Pyruvate A (from C1-C2-C3) -> Its carboxyl carbon comes from C1 of glucose.")
    print("Pyruvate B (from C4-C5-C6) -> Its carboxyl carbon comes from C4 of glucose.\n")

    # Step 3: Simulate pyruvate decarboxylation to release CO2.
    # The carboxyl carbon (C1) of each pyruvate is released as CO2.
    print("--- Step 3: Pyruvate Decarboxylation ---")
    print("Each pyruvate molecule releases its carboxyl carbon as CO2.")

    labeled_co2_count = 0
    num_from_A = 0
    num_from_B = 0
    
    # Check Pyruvate A
    if pyruvate_A_carboxyl_label == 13:
        labeled_co2_count += 1
        num_from_A = 1
        print("Result for Pyruvate A: Carboxyl group is 13C. One molecule of 13CO2 is released.")
    else:
        print("Result for Pyruvate A: Carboxyl group is 12C. One molecule of 12CO2 is released.")

    # Check Pyruvate B
    if pyruvate_B_carboxyl_label == 13:
        labeled_co2_count += 1
        num_from_B = 1
        print("Result for Pyruvate B: Carboxyl group is 13C. One molecule of 13CO2 is released.\n")
    else:
        print("Result for Pyruvate B: Carboxyl group is 12C. One molecule of 12CO2 is released.\n")

    # Step 4: Final calculation and summary.
    print("--- Step 4: Final Calculation ---")
    print(f"The final count of 13CO2 molecules is the sum of those released from each pyruvate.")
    print(f"Equation: {num_from_A} (from Pyruvate A) + {num_from_B} (from Pyruvate B) = {labeled_co2_count}")
    print(f"\nTherefore, a total of {labeled_co2_count} 13C-labeled CO2 molecules will be released.")

trace_glycolysis_co2_release()