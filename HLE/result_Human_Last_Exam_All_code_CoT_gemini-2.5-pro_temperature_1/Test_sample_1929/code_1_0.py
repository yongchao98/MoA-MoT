def identify_reaction_product():
    """
    This script identifies the product of the reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene.
    """
    # Define reactants and product formula
    diene = "Butadiene (C4H6)"
    dienophile = "1,1-dichloro-2,2-difluoroethene (C2Cl2F2)"
    product_formula = "C6H6Cl2F2"
    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"

    # Explain the reaction
    print("1. The reaction is a Diels-Alder [4+2] cycloaddition.")
    print(f"   - Diene: {diene}")
    print(f"   - Dienophile: {dienophile}")

    print("\n2. A six-membered ring is formed from the 4 carbons of the diene and the 2 carbons of the dienophile.")
    
    # Show the overall reaction equation
    print("\n3. The final chemical equation is:")
    print(f"   {diene} + {dienophile} -> {product_name} ({product_formula})")

    # State the product name clearly
    print(f"\n4. The name of the product is: {product_name}")

    # As requested, output each number from the product's name (the final "equation")
    print("\n5. The numbers in the IUPAC name of the product are:")
    numbers_in_name = [4, 4, 5, 5, 1]
    for number in numbers_in_name:
        print(number)

if __name__ == "__main__":
    identify_reaction_product()