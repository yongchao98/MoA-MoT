def solve_glycolysis_labeling():
    """
    Traces 13C labels from 1,4-13C glucose through glycolysis to determine
    the number of labeled CO2 molecules produced.
    """
    
    # --- Step 1 & 2: Trace the labels from Glucose to Pyruvate ---
    print("This script traces the path of labeled carbons from 1,4-13C glucose through glycolysis.")
    print("-" * 70)

    # Initial state: Glucose is a 6-carbon molecule.
    # We represent it as a list where the index corresponds to the carbon number (1-based).
    # '13C' indicates a labeled carbon, '12C' is a standard carbon.
    glucose = ['13C', '12C', '12C', '13C', '12C', '12C']
    print(f"Starting molecule: 1,4-13C glucose. Carbons [C1, C4] are labeled.")
    print(f"Glucose: C1({glucose[0]})-C2({glucose[1]})-C3({glucose[2]})-C4({glucose[3]})-C5({glucose[4]})-C6({glucose[5]})")
    print("\nGlycolysis begins. The 6-carbon glucose is cleaved into two 3-carbon molecules.")
    
    # Fructose-1,6-bisphosphate cleaves into DHAP (from C1, C2, C3) and G3P (from C4, C5, C6)
    dhap = [glucose[0], glucose[1], glucose[2]] # Carbons from glucose C1, C2, C3
    g3p_from_glucose_bottom = [glucose[3], glucose[4], glucose[5]] # Carbons from glucose C4, C5, C6
    
    print("-> Cleavage results in:")
    print(f"   1. Dihydroxyacetone Phosphate (DHAP) from C1-C3 -> Labeled at its C1 position (from Glucose C1)")
    print(f"   2. Glyceraldehyde-3-Phosphate (G3P) from C4-C6 -> Labeled at its C1 position (from Glucose C4)")
    
    # DHAP is isomerized to G3P. This isomerization inverts the carbon chain.
    # DHAP C1 becomes G3P C3. DHAP C3 becomes G3P C1.
    g3p_from_dhap = [dhap[2], dhap[1], dhap[0]]

    print("\n-> DHAP is then isomerized to form a second molecule of G3P.")
    print(f"   The G3P derived from DHAP is now labeled at its C3 position.")

    # We now have two G3P molecules ready for the payoff phase.
    # g3p_1 is from the bottom half of glucose. Its label is at C1.
    # g3p_2 is from the top half (via DHAP). Its label is at C3.
    
    # Convert G3P to Pyruvate.
    # G3P C1 (aldehyde) becomes Pyruvate C3 (methyl).
    # G3P C3 (phosphate) becomes Pyruvate C1 (carboxyl).
    pyruvate_1_label = g3p_from_glucose_bottom[0] # The label from G3P C1 ends up on Pyruvate C3.
    pyruvate_2_label = g3p_from_dhap[2]         # The label from G3P C3 ends up on Pyruvate C1.

    print("\n-> The two G3P molecules are converted to two Pyruvate molecules:")
    print(f"   - The first G3P produces Pyruvate labeled at the C3 (methyl) position.")
    print(f"   - The second G3P produces Pyruvate labeled at the C1 (carboxyl) position.")

    # --- Step 3 & 4: Analyze CO2 Production ---
    print("\n" + "-" * 70)
    print("Step 3: Analyze CO2 Production from Glycolysis")
    print("-" * 70)

    print("The metabolic pathway of glycolysis is defined as the sequence of reactions that converts glucose into two molecules of pyruvate.")
    print("In this entire 10-step process, the carbon backbone is converted from one 6-carbon molecule into two 3-carbon molecules.")
    print("\nCrucially, NO carbon atoms are lost or released as carbon dioxide (CO2) during glycolysis.")

    num_labeled_co2 = 0
    print(f"Therefore, the number of 13C-labeled CO2 molecules released is: {num_labeled_co2}")

    print("\n(Note: CO2 is released in subsequent pathways, like the Pyruvate Dehydrogenase Complex reaction, but not in glycolysis itself.)")

    # --- Step 5: Final Equation ---
    print("\n" + "-" * 70)
    print("Final Equation")
    print("-" * 70)
    
    num_glucose = 1
    num_pyruvate = 2
    
    print(f"The reaction for glycolysis shows:")
    print(f"{num_glucose} Glucose -> {num_pyruvate} Pyruvate + {num_labeled_co2} CO2")
    
solve_glycolysis_labeling()