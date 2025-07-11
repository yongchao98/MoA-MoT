def find_elimination_product():
    """
    Analyzes the E2 elimination of (1S,2R)-1-bromo-2-methylcyclohexane
    and determines the name of the product.
    """
    print("Task: Identify the product of (1S,2R)-1-bromo-2-methylcyclohexane treated with potassium tert-butoxide.")
    print("-" * 70)
    
    print("Step 1: Analyzing the reaction conditions.")
    print("-> Substrate: (1S,2R)-1-bromo-2-methylcyclohexane (a secondary alkyl halide).")
    print("-> Reagent: Potassium tert-butoxide (a strong, sterically hindered base).")
    print("-> Conclusion: The reaction is an E2 (bimolecular elimination) reaction.")
    print("-" * 70)

    print("Step 2: Understanding the stereochemical requirements for an E2 reaction.")
    print("-> The E2 mechanism requires an 'anti-periplanar' geometry.")
    print("-> In a cyclohexane chair, this means the leaving group (Br) and the beta-proton (H) MUST BOTH BE AXIAL.")
    print("-" * 70)

    print("Step 3: Analyzing the conformations of the starting material.")
    print("-> (1S,2R)-1-bromo-2-methylcyclohexane is a trans-1,2-disubstituted cyclohexane.")
    print("-> It exists primarily in a stable 'diequatorial' conformation (Br and methyl are both equatorial).")
    print("-> This stable conformer is UNREACTIVE for E2 because the Br is not axial.")
    print("-> The reaction must proceed through the less stable 'diaxial' conformer, where Br and methyl are both axial.")
    print("-" * 70)

    print("Step 4: Identifying the site of elimination from the reactive (diaxial) conformer.")
    print("-> In the diaxial conformer, the Br is AXIAL on C1.")
    print("-> We look for an AXIAL proton on an adjacent carbon (C2 or C6).")
    print("  - On C2: The methyl group is AXIAL, so the proton on C2 is EQUATORIAL. Elimination here is IMPOSSIBLE.")
    print("  - On C6: This carbon has an AXIAL proton which is anti-periplanar to the axial Br.")
    print("-> Therefore, the base must abstract the axial proton from C6.")
    print("-" * 70)
    
    print("Step 5: Naming the final product.")
    print("-> The elimination forms a double bond between C1 and C6.")
    print("-> The methyl group on C2 remains.")
    print("-> Numbering the new alkene ring gives the methyl group at position 3.")
    print("The final product is named based on these parts:")
    
    product_parts = {
        "Position": "3",
        "Separator": "-",
        "Substituent": "methyl",
        "Ring System": "cyclohexene"
    }

    final_name = ""
    for key, value in product_parts.items():
        print(f"  - {key}: {value}")
        final_name += value
    
    print("\n------------------- FINAL PRODUCT -------------------")
    print(f"The name of the major product is: {final_name}")
    print("-----------------------------------------------------")

find_elimination_product()