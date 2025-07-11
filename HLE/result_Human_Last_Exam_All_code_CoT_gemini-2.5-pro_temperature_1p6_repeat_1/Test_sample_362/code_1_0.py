def solve_wittig_reaction():
    """
    Determines the product of a specified Wittig reaction and prints the details.
    """
    # 1. Define Reactants
    aldehyde = {
        "name": "Pivalaldehyde (2,2-dimethylpropanal)",
        "structure": "(CH3)3C-CH=O",
        "formula": "C5H10O"
    }

    wittig_reagent = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "ylide_part": "=CH-CH2-(C6H4Cl)" # Where C6H4Cl is a 2-chlorophenyl group
    }

    # 2. Define Products by applying the Wittig reaction mechanism
    # The (CH3)3C-CH= part from the aldehyde combines with the =CH-CH2-(C6H4Cl) part from the ylide.
    main_product = {
        "name": "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene",
        "structure": "(CH3)3C-CH=CH-CH2-(C6H4Cl)",
        "formula": "C13H17Cl" # C(4+1+1+1+6) H(9+1+1+2+4) Cl(1)
    }

    byproduct = {
        "name": "Triphenylphosphine oxide",
        "structure": "(C6H5)3P=O",
        "formula": "C18H15OP"
    }

    # 3. Print the reaction explanation and results
    print("--- Wittig Reaction Analysis ---")
    print("\nStep 1: Identify the reactants.")
    print(f"Aldehyde: {aldehyde['name']}")
    print(f"Structure Fragment: {aldehyde['structure']}")
    print(f"Wittig Reagent: {wittig_reagent['name']}")
    print(f"Ylide Fragment: Ph3P{wittig_reagent['ylide_part']}")

    print("\nStep 2: Combine the fragments to form the product.")
    print("The C=O bond of the aldehyde and the P=C bond of the ylide are broken.")
    print("A new C=C bond is formed, along with a P=O byproduct.")

    print("\n--- Final Reaction Equation ---")
    # To "output each number", we display the molecular formulas which contain the atom counts.
    print(f"{aldehyde['structure']}   +   Ph3P{wittig_reagent['ylide_part']}   --->   {main_product['structure']}   +   {byproduct['structure']}")
    print(f"({aldehyde['formula']})        ({wittig_reagent['name']})           ({main_product['formula']})          ({byproduct['formula']})")

    print("\n--- Product Details ---")
    print(f"The main product is: {main_product['name']}")
    print(f"Molecular Formula: {main_product['formula']}")

if __name__ == '__main__':
    solve_wittig_reaction()