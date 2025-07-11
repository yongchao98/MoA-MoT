import re

def chemical_reaction_product():
    """
    Analyzes the Diels-Alder reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the details.
    """
    # Define reactants and product information
    diene = {
        "name": "Butadiene",
        "formula": "C4H6"
    }
    dienophile = {
        "name": "1,1-dichloro-2,2-difluoroethene",
        "formula": "C2Cl2F2"
    }
    product = {
        "name": "4,4-dichloro-5,5-difluorocyclohexene",
        "formula": "C6H6Cl2F2",
        "smiles": "C1=CCC(Cl)(Cl)C(F)(F)C1"
    }

    print("### Analysis of the Chemical Reaction ###")
    print("\nThe reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")
    print("It forms a six-membered ring, resulting in a cyclohexene derivative.")

    # Print the reaction equation
    print("\nReaction Equation:")
    final_equation = f"{diene['formula']} + {dienophile['formula']} -> {product['formula']}"
    print(final_equation)
    print(f"{diene['name']} + {dienophile['name']} -> {product['name']}")

    # Print product details
    print("\nProduct Details:")
    print(f"  IUPAC Name: {product['name']}")
    print(f"  Chemical Formula: {product['formula']}")
    print(f"  SMILES String: {product['smiles']}")

    # Fulfill the special request to output numbers from the final equation
    print("\nNumbers from the final equation's chemical formulas:")
    print(f"Equation: {final_equation}")

    def extract_and_print_numbers(label, formula):
        """Helper function to extract and print numbers from a formula string."""
        # A '1' is implied for elements with no subscript. We will only extract explicit numbers.
        numbers = re.findall(r'\d+', formula)
        print(f"Numbers in {label} ({formula}): ", end="")
        for num in numbers:
            print(num, end=" ")
        print()

    extract_and_print_numbers("Reactant 1", diene['formula'])
    extract_and_print_numbers("Reactant 2", dienophile['formula'])
    extract_and_print_numbers("Product", product['formula'])

# Execute the function
chemical_reaction_product()