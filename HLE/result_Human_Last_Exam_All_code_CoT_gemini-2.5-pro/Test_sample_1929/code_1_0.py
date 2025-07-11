def solve_reaction():
    """
    Determines and describes the product of the reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    # Define Reactants
    reactant1_name = "Butadiene"
    reactant1_formula_parts = {"C": 4, "H": 6}
    
    reactant2_name = "1,1-dichloro-2,2-difluoroethene"
    reactant2_formula_parts = {"C": 2, "Cl": 2, "F": 2}

    # Define Product based on Diels-Alder cycloaddition
    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    product_formula_parts = {
        "C": reactant1_formula_parts["C"] + reactant2_formula_parts["C"],
        "H": reactant1_formula_parts["H"],
        "Cl": reactant2_formula_parts["Cl"],
        "F": reactant2_formula_parts["F"]
    }

    # Print the analysis
    print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder [4+2] cycloaddition.")
    print("-" * 80)
    
    print(f"Reactant 1 (Diene): {reactant1_name}")
    print(f"Reactant 2 (Dienophile): {reactant2_name}")
    print("\nIn this reaction, the diene and dienophile combine to form a new six-membered ring.")
    
    print(f"\nThe chemical product is: {product_name}")
    print("-" * 80)
    
    # Print the final chemical equation, outputting each number as requested
    print("The final balanced chemical equation is:")
    
    # C4H6 + C2Cl2F2 -> C6H6Cl2F2
    print(
        f"C{reactant1_formula_parts['C']}H{reactant1_formula_parts['H']} + "
        f"C{reactant2_formula_parts['C']}Cl{reactant2_formula_parts['Cl']}F{reactant2_formula_parts['F']} -> "
        f"C{product_formula_parts['C']}H{product_formula_parts['H']}Cl{product_formula_parts['Cl']}F{product_formula_parts['F']}"
    )

solve_reaction()