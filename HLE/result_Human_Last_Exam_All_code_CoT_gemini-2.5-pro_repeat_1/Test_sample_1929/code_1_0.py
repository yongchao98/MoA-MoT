def solve_reaction():
    """
    Identifies the product of the reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and displays the reaction details.
    """
    # 1. Define reactants and reaction type
    reactant1_name = "butadiene"
    reactant1_formula = "C4H6"

    reactant2_name = "1,1-dichloro-2,2-difluoroethene"
    reactant2_formula = "C2Cl2F2"

    reaction_type = "Diels-Alder [4+2] cycloaddition"

    # 2. Determine the product
    # The reaction forms a six-membered ring called cyclohexene.
    # The substituents from the dienophile remain on the carbons that were part of its double bond.
    # The resulting product is named 4,4-dichloro-5,5-difluorocyclohexene.
    product_name = "4,4-dichloro-5,5-difluorocyclohexene"
    # The product formula is the sum of the reactant atoms: C(4+2)H6Cl2F2
    product_formula = "C6H6Cl2F2"

    # 3. Define the balanced equation's coefficients
    coeff_reactant1 = 1
    coeff_reactant2 = 1
    coeff_product = 1

    # 4. Print the detailed explanation and results
    print(f"The reaction between {reactant1_name} and {reactant2_name} is a {reaction_type}.")
    print("A new six-membered ring is formed, resulting in a cyclohexene derivative.")
    print("\n--- Reaction Summary ---")
    print(f"Product Name: {product_name}")
    print(f"Product Formula: {product_formula}")

    print("\n--- Balanced Chemical Equation ---")
    # Print the full equation string
    print(f"{coeff_reactant1} {reactant1_formula} + {coeff_reactant2} {reactant2_formula} -> {coeff_product} {product_formula}")

    print("\nThe numbers (stoichiometric coefficients) in the final equation are:")
    # As requested, output each number from the equation
    print(coeff_reactant1)
    print(coeff_reactant2)
    print(coeff_product)

# Execute the function
solve_reaction()