def solve_reaction():
    """
    Identifies the product of the reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the balanced chemical equation.
    """

    # Define the reactants and their properties
    diene = {
        "name": "1,3-butadiene",
        "formula": "C4H6",
        "coeff": 1
    }

    dienophile = {
        "name": "1,1-dichloro-2,2-difluoroethene",
        "formula": "C2Cl2F2",
        "coeff": 1
    }

    # Determine the product of the Diels-Alder reaction
    product = {
        "name": "3,3-dichloro-4,4-difluorocyclohexene",
        "formula": "C6H6Cl2F2",
        "coeff": 1
    }

    # Print the details of the reaction
    print("Reaction Type: Diels-Alder ([4+2] Cycloaddition)")
    print(f"Product Name: {product['name']}")
    print("\nBalanced Chemical Equation:")

    # Print the final equation, including each number (stoichiometric coefficient)
    print(f"{diene['coeff']} {diene['formula']} + {dienophile['coeff']} {dienophile['formula']} -> {product['coeff']} {product['formula']}")

# Execute the function
solve_reaction()