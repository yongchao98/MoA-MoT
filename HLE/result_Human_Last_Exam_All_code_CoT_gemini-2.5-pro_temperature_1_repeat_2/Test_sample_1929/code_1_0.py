def solve_reaction():
    """
    This script identifies the product of the reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the details.
    """

    # Define reactant and product information
    diene_name = "Butadiene"
    diene_formula = "C4H6"
    dienophile_name = "1,1-dichloro-2,2-difluoroethene"
    dienophile_formula = "C2Cl2F2"
    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    product_formula = "C6H6Cl2F2"

    # Explain the reaction
    print(f"The reaction between {diene_name} and {dienophile_name} is a [4+2] cycloaddition (Diels-Alder reaction).")
    print("This reaction forms a six-membered ring.")
    print("-" * 60)

    # Print the reaction equation, which includes the required numbers (subscripts)
    print("Final Reaction Equation:")
    # The numbers in the equation are the subscripts for each element.
    # For C4H6: 4, 6
    # For C2Cl2F2: 2, 2, 2
    # For C6H6Cl2F2: 6, 6, 2, 2
    print(f"    {diene_formula} + {dienophile_formula} -> {product_formula}")
    print("-" * 60)

    # Print the product name, which includes the required numbers (locants)
    print("The IUPAC name of the product is:")
    # The numbers in the name are the locants for the substituent groups and the double bond.
    # The locants are 4, 4, 5, 5, and 1.
    print(f"    {product_name}")

if __name__ == "__main__":
    solve_reaction()