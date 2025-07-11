def trace_glucose_catabolism():
    """
    Traces the carbons of 1,4-13C glucose through glycolysis and pyruvate
    decarboxylation to determine how many labeled CO2 molecules are released.
    """
    print("--- Tracing 1,4-13C Glucose Catabolism ---")

    # Step 1: Define the starting labeled glucose molecule
    # C1 and C4 are labeled with 13C.
    # We use a dictionary to map carbon position to its label status.
    glucose = {1: '13C', 2: 'C', 3: 'C', 4: '13C', 5: 'C', 6: 'C'}
    print(f"\nInitial Molecule: 1,4-13C Glucose")
    print(f"Carbon labels: {glucose}")
    print("-" * 40)

    # Step 2: Glycolysis
    # Glucose (6C) is cleaved into DHAP (from C1,C2,C3) and G3P (from C4,C5,C6).
    print("Step A: Glycolysis - Cleavage")
    dhap = { 'C1': glucose[1], 'C2': glucose[2], 'C3': glucose[3] }
    g3p_1 = { 'C1': glucose[4], 'C2': glucose[5], 'C3': glucose[6] }
    print(f"Glucose is split into:")
    print(f"  - DHAP (from Glucose C1,C2,C3) -> {dhap}")
    print(f"  - G3P (from Glucose C4,C5,C6)  -> {g3p_1}")

    # DHAP is isomerized to G3P. Now we have two G3P molecules.
    g3p_2 = dhap
    print("\nDHAP is isomerized to G3P, resulting in two G3P molecules.")
    print(f"  - G3P_1 (original) -> {g3p_1}")
    print(f"  - G3P_2 (from DHAP) -> {g3p_2}")
    print("-" * 40)

    # Both G3P molecules are converted to Pyruvate.
    # Carbon Mapping: G3P C1 -> Pyruvate C3 (carboxyl), G3P C2 -> Pyruvate C2, G3P C3 -> Pyruvate C1
    print("Step B: Glycolysis - Payoff Phase")
    pyruvate_1 = {'methyl_C': g3p_1['C3'], 'keto_C': g3p_1['C2'], 'carboxyl_C': g3p_1['C1']}
    pyruvate_2 = {'methyl_C': g3p_2['C3'], 'keto_C': g3p_2['C2'], 'carboxyl_C': g3p_2['C1']}
    print("The two G3P molecules are converted to two Pyruvate molecules:")
    print(f"  - Pyruvate_1 (from G3P_1): {pyruvate_1}")
    print(f"  - Pyruvate_2 (from G3P_2): {pyruvate_2}")
    print("-" * 40)

    # Step 3: Pyruvate Decarboxylation
    # The carboxyl_C of each pyruvate is released as CO2.
    print("Step C: Pyruvate Decarboxylation (Link Reaction)")
    co2_from_pyruvate1 = pyruvate_1['carboxyl_C']
    co2_from_pyruvate2 = pyruvate_2['carboxyl_C']
    print("The carboxyl carbon from each pyruvate is released as CO2.")
    print(f"  - CO2 from Pyruvate_1 originates from Glucose C4 -> '{co2_from_pyruvate1}'")
    print(f"  - CO2 from Pyruvate_2 originates from Glucose C1 -> '{co2_from_pyruvate2}'")
    print("-" * 40)

    # Step 4: Count the labeled CO2 molecules
    print("Final Calculation:")
    labeled_co2_count = 0
    calculation_terms = []

    if '13C' in co2_from_pyruvate1:
        labeled_co2_count += 1
        calculation_terms.append("1")
    else:
        calculation_terms.append("0")

    if '13C' in co2_from_pyruvate2:
        labeled_co2_count += 1
        calculation_terms.append("1")
    else:
        calculation_terms.append("0")

    print("We sum the number of labeled CO2 molecules produced:")
    equation = f"{calculation_terms[0]} + {calculation_terms[1]} = {labeled_co2_count}"
    print(f"Equation: {equation}")
    print(f"\nConclusion: A total of {labeled_co2_count} molecule(s) of 13C-labeled CO2 are released.")

if __name__ == '__main__':
    trace_glucose_catabolism()

<<<2>>>