def solve_wittig_reaction():
    """
    Determines and displays the product of a Wittig reaction between
    pivalaldehyde and (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane.
    """

    # 1. Define Reactants
    pivalaldehyde = {
        "name": "Pivalaldehyde",
        "formula": "C5H10O",
        "structure": "(CH3)3CCHO",
        "smiles": "CC(C)(C)C=O"
    }

    wittig_reagent = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "formula": "C26H22ClP",
        "structure": "Ph3P=CHCH2(C6H4Cl)",
        "smiles": "Clc1ccccc1CC=[P](c2ccccc2)(c3ccccc3)c4ccccc4"
    }

    # 2. Determine Products based on the Wittig reaction mechanism
    # The C=O group of the aldehyde is replaced by the C=CH-R group from the ylide.
    # Product: (CH3)3C-CH=CH-CH2-(C6H4Cl)
    # Byproduct: Ph3P=O
    main_product = {
        "name": "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene",
        "formula": "C13H17Cl",
        "structure": "(CH3)3CCH=CHCH2(C6H4Cl)",
        "smiles": "CC(C)(C)C=CCc1c(Cl)cccc1"
    }

    byproduct = {
        "name": "Triphenylphosphine oxide",
        "formula": "C18H15OP",
        "structure": "Ph3PO",
        "smiles": "O=P(c1ccccc1)(c2ccccc2)c3ccccc3"
    }
    
    # 3. Print the reaction equation
    print("The Wittig Reaction:")
    print("This reaction forms a C=C double bond by replacing a carbonyl C=O bond with a carbon-phosphorus ylide C=P bond.")
    print("-" * 50)
    print("Reactants:")
    print(f"  Aldehyde: {pivalaldehyde['name']}")
    print(f"    - Formula: {pivalaldehyde['formula']}")
    print(f"    - SMILES: {pivalaldehyde['smiles']}")
    print(f"  Wittig Reagent: {wittig_reagent['name']}")
    print(f"    - Formula: {wittig_reagent['formula']}")
    print(f"    - SMILES: {wittig_reagent['smiles']}")
    print("-" * 50)
    
    print("The Final Equation:")
    # Using 'structure' for a simplified view of the equation
    # The numbers in the formulas are part of the output as requested.
    reactant1_str = f"{pivalaldehyde['name']} ({pivalaldehyde['formula']})"
    reactant2_str = f"{wittig_reagent['name']} ({wittig_reagent['formula']})"
    product_str = f"{main_product['name']} ({main_product['formula']})"
    byproduct_str = f"{byproduct['name']} ({byproduct['formula']})"
    
    print(f"{reactant1_str} + {reactant2_str} -> {product_str} + {byproduct_str}")
    
    print("-" * 50)
    print("Product Details:")
    print(f"  Main Product: {main_product['name']}")
    print(f"    - Formula: {main_product['formula']}")
    print(f"    - SMILES: {main_product['smiles']}")
    print(f"  Byproduct: {byproduct['name']}")
    print(f"    - Formula: {byproduct['formula']}")
    print(f"    - SMILES: {byproduct['smiles']}")


solve_wittig_reaction()