def solve_reaction():
    """
    Identifies the product of the reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the details.
    """
    # Define reactant and product information
    reactant1 = {"name": "Butadiene", "formula": "C4H6"}
    reactant2 = {"name": "1,1-dichloro-2,2-difluoroethene", "formula": "C2Cl2F2"}
    product = {"name": "4,4-dichloro-5,5-difluorocyclohexene", "formula": "C6H6Cl2F2"}
    reaction_type = "Diels-Alder [4+2] cycloaddition"

    # Explain the reaction
    print(f"The reaction between {reactant1['name']} and {reactant2['name']} is a {reaction_type}.")
    print("This reaction forms a six-membered ring, resulting in a single product.")
    print("-" * 60)

    # Print the final reaction equation
    print("The final reaction equation is:")
    print(f"{reactant1['name']} + {reactant2['name']} -> {product['name']}")
    print(f"{reactant1['formula']} + {reactant2['formula']} -> {product['formula']}")
    print("-" * 60)

    # Fulfill the "output each number in the final equation" requirement
    # The stoichiometric coefficients are all 1. The numbers are from the molecular formulas.
    print("The numbers in the final chemical equation (from the molecular formulas) are:")
    print(f"Reactant 1 ({reactant1['formula']}): 4, 6")
    print(f"Reactant 2 ({reactant2['formula']}): 2, 2, 2")
    print(f"Product ({product['formula']}): 6, 6, 2, 2")

if __name__ == "__main__":
    solve_reaction()