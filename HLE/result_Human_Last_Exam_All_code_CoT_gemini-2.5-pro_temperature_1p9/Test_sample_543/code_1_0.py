def find_product_name():
    """
    This script determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """

    # Reactants and conditions
    reactant1 = "1,3-dibromo-2-iodobenzene"
    reactant2 = "excess phenyl magnesium bromide (PhMgBr)"
    conditions = "refluxing in THF, followed by aqueous work-up"

    print(f"Starting reaction of {reactant1} with {reactant2} under {conditions}.\n")

    print("Step 1: The most reactive halogen, Iodine, is substituted first.")
    print("Reaction: 1,3-dibromo-2-iodobenzene + PhMgBr -> 1,3-dibromo-2-phenylbenzene\n")

    print("Step 2: Due to excess Grignard reagent and reflux, the less reactive Bromine atoms are also substituted.")
    print("Reaction: 1,3-dibromo-2-phenylbenzene + 2 PhMgBr -> 1,2,3-triphenylbenzene\n")

    final_product = {
        "name": "1,2,3-triphenylbenzene",
        "numbers": [1, 2, 3]
    }

    print("The final organic product has three phenyl groups on a central benzene ring.")
    print("The IUPAC name is based on the positions of these phenyl groups.")
    
    # Building and printing the final name from its components as requested
    name_parts = [str(n) for n in final_product["numbers"]]
    final_name_string = ",".join(name_parts) + "-triphenylbenzene"
    
    print(f"\nThe IUPAC name of the final product is: {final_name_string}")

if __name__ == "__main__":
    find_product_name()