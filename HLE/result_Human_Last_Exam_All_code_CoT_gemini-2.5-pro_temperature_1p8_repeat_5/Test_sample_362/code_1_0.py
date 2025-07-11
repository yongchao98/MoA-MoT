import re

def solve_wittig_reaction():
    """
    This function determines and prints the product of the specified Wittig reaction.
    """
    # 1. Define the reactants
    aldehyde = "pivalaldehyde [(CH3)3C-CHO]"
    ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane [(2-Cl-C6H4)-CH2-CH=PPh3]"

    # 2. Determine the product's structure and name
    # The C=O of the aldehyde and the C=P of the ylide form a new C=C bond.
    product_structure = "(CH3)3C-CH=CH-CH2-(C6H4Cl)"
    product_byproduct = "triphenylphosphine oxide (O=PPh3)"
    
    # Unstabilized ylides typically give the Z-alkene as the major product.
    product_iupac_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 3. Print the results
    print("The Wittig reaction occurs between:")
    print(f"Reactant 1 (Aldehyde): {aldehyde}")
    print(f"Reactant 2 (Ylide): {ylide}")
    print("-" * 20)
    print(f"The reaction forms a new alkene and a byproduct.")
    print(f"\nProduct Structure: {product_structure}")
    print(f"Byproduct: {product_byproduct}")
    
    print("\nThe IUPAC name for the major product is:")
    print(product_iupac_name)

    # 4. As requested, output each number from the final product name
    numbers = re.findall(r'\d+', product_iupac_name)
    print("\nThe numbers found in the product's IUPAC name are:")
    for num in numbers:
        print(num)

if __name__ == "__main__":
    solve_wittig_reaction()