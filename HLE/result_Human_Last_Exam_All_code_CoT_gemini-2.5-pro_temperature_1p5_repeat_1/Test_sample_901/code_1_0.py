import textwrap

def solve_reaction():
    """
    Analyzes the E2 elimination of (1S,2R)-1-bromo-2-methylcyclohexane
    and identifies the major product.
    """
    
    # Step 1: Define the reactants and reaction conditions.
    starting_material = "(1S,2R)-1-bromo-2-methylcyclohexane"
    reagent = "Potassium tert-butoxide (a strong, bulky base)"
    reaction_type = "E2 Elimination"

    # Step 2: Explain the mechanism and stereochemical requirements.
    explanation = f"""
    The reaction of {starting_material} with {reagent} is an {reaction_type}.
    This reaction mechanism requires an anti-periplanar geometry, meaning the leaving group (Br) and a beta-hydrogen (H on an adjacent carbon) must both be in axial positions on the cyclohexane chair.
    """

    # Step 3: Analyze the chair conformations of the starting material.
    explanation += """
    The starting material, (1S,2R)-1-bromo-2-methylcyclohexane, has the bromo and methyl groups cis to each other. This leads to two possible chair conformations in equilibrium:
    Conformer A: Axial Bromo, Equatorial Methyl
    Conformer B: Equatorial Bromo, Axial Methyl

    The methyl group is bulkier than the bromo group. To minimize steric strain, the bulkier group prefers the equatorial position. Therefore, Conformer A is the more stable and predominant conformer. The reaction proceeds primarily through Conformer A.
    """

    # Step 4: Identify possible products based on the reactive conformer.
    explanation += """
    In Conformer A, the bromo group is axial, which is required for the E2 reaction. We now look for axial beta-hydrogens on the adjacent carbons (C2 and C6).

    1. Elimination at C2: The H on C2 is axial. Removing this H gives the Zaitsev product: 1-methylcyclohexene.
       (Product 1 Name: 1-methylcyclohexene)

    2. Elimination at C6: The H on C6 is also axial. Removing this H gives the Hofmann product: 3-methylcyclohexene.
       (Product 2 Name: 3-methylcyclohexene)
    """

    # Step 5: Determine the major product based on the bulky base.
    explanation += """
    Potassium tert-butoxide is a sterically hindered (bulky) base. Bulky bases preferentially remove the less sterically hindered proton. The proton on C6 is more accessible than the proton on C2 (which is shielded by the methyl group).

    Therefore, the reaction favors the formation of the less substituted alkene, the Hofmann product.
    """
    
    # Final conclusion
    final_product_name = "3-methylcyclohexene"
    conclusion = f"\nThe major product of the reaction is {final_product_name}."

    # Print the step-by-step explanation and the final answer.
    print(textwrap.dedent(explanation))
    print("="*40)
    print("CONCLUSION")
    print("="*40)
    print(conclusion)

# Execute the function to find the answer
solve_reaction()