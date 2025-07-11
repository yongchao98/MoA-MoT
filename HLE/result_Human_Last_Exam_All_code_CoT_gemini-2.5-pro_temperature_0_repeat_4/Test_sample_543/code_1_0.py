def get_product_name():
    """
    Analyzes the chemical reaction and provides the IUPAC name of the product.
    """
    # 1. Define reactants and conditions
    substrate = "1,3-dibromo-2-iodobenzene"
    reagent = "excess phenyl magnesium bromide (PhMgBr)"
    
    print(f"Analyzing the reaction of {substrate} with {reagent}.")
    print("-" * 60)

    # 2. Explain the chemical principles
    print("Principle of Reactivity:")
    print("The reaction involves a Grignard reagent substituting halogens on a benzene ring.")
    print("The reactivity of the carbon-halogen bonds follows the order: C-I > C-Br > C-Cl.")
    print("This means the iodine atom will react first, followed by the bromine atoms.")
    print("-" * 60)

    # 3. Describe the reaction sequence
    print("Reaction Sequence:")
    print("a) First, the highly reactive phenyl magnesium bromide replaces the most labile halogen, the iodine atom at position 2, with a phenyl group.")
    print("b) Because an *excess* of the Grignard reagent is used and the reaction is heated (reflux), the reaction is driven further.")
    print("c) The two less reactive bromine atoms at positions 1 and 3 are also subsequently replaced by phenyl groups.")
    print("-" * 60)

    # 4. Identify the final product and its IUPAC name
    print("Final Product and IUPAC Name:")
    print("The final product is a benzene ring where the substituents at positions 1, 2, and 3 are all phenyl groups.")
    
    final_product_name = "1,2,3-triphenylbenzene"
    
    print(f"\nThe IUPAC name of this product is: {final_product_name}")
    
    # As requested, explicitly outputting the numbers involved in the name
    print("\nThe numbers in the name, 1, 2, and 3, specify the positions of the three 'phenyl' groups on the central 'benzene' ring.")
    print("-" * 60)

if __name__ == '__main__':
    get_product_name()