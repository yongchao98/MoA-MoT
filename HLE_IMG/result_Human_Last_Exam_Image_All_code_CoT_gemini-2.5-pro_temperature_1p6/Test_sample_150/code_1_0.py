def track_reaction_sequence():
    """
    This function tracks a multi-step organic synthesis reaction,
    identifying each intermediate and the final product.
    """
    print("Tracking the chemical reaction sequence step-by-step:\n")

    # Step 1: Friedel-Crafts Acylation
    step1 = {
        "name": "Step 1: Friedel-Crafts Acylation",
        "reactants": "Benzene and Propanoyl chloride",
        "reagents": "AlCl3",
        "product_name": "Propiophenone (1-phenylpropan-1-one)",
        "explanation": "The propanoyl group attaches to the benzene ring to form a ketone."
    }
    print(f"{step1['name']}")
    print(f"  Reactants: {step1['reactants']}")
    print(f"  Reagents: {step1['reagents']}")
    print(f"  Explanation: {step1['explanation']}")
    print(f"  Product (Intermediate-1): {step1['product_name']}\n")
    intermediate_1 = step1['product_name']

    # Step 2: Electrophilic Aromatic Bromination
    step2 = {
        "name": "Step 2: Electrophilic Aromatic Bromination",
        "reactant": intermediate_1,
        "reagents": "Br2/FeBr3",
        "product_name": "1-(3-bromophenyl)propan-1-one",
        "explanation": "The acyl group is a meta-director, so bromine adds to the meta (position 3) on the phenyl ring."
    }
    print(f"{step2['name']}")
    print(f"  Reactant: {step2['reactant']}")
    print(f"  Reagents: {step2['reagents']}")
    print(f"  Explanation: {step2['explanation']}")
    print(f"  Product (Intermediate-2): {step2['product_name']}\n")
    intermediate_2 = step2['product_name']

    # Step 3: Catalytic Hydrogenation (Hydrogenolysis)
    step3 = {
        "name": "Step 3: Catalytic Hydrogenation",
        "reactant": intermediate_2,
        "reagents": "H2/Pd",
        "product_name": "1-(3-bromophenyl)propane",
        "explanation": "The benzylic ketone C=O group is completely reduced to a CH2 group."
    }
    print(f"{step3['name']}")
    print(f"  Reactant: {step3['reactant']}")
    print(f"  Reagents: {step3['reagents']}")
    print(f"  Explanation: {step3['explanation']}")
    print(f"  Product (Intermediate-3): {step3['product_name']}\n")
    intermediate_3 = step3['product_name']

    # Step 4: Radical Benzylic Bromination
    step4 = {
        "name": "Step 4: Radical Benzylic Bromination",
        "reactant": intermediate_3,
        "reagents": "NBS, (PhCO2)2, CCl4",
        "product_name": "1-bromo-1-(3-bromophenyl)propane",
        "explanation": "NBS selectively brominates the benzylic position (carbon 1 of the propane chain)."
    }
    print(f"{step4['name']}")
    print(f"  Reactant: {step4['reactant']}")
    print(f"  Reagents: {step4['reagents']}")
    print(f"  Explanation: {step4['explanation']}")
    print(f"  Final Product: {step4['product_name']}\n")
    final_product = step4['product_name']

    print("---CONCLUSION---")
    print(f"The name of the final product is: {final_product}")

if __name__ == "__main__":
    track_reaction_sequence()