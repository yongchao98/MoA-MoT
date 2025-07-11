import sys

def solve_reaction():
    """
    This function determines the IUPAC name of the product from the specified reaction.
    """
    # --- Reaction Analysis ---
    reactant = "1,3-dibromo-2-iodobenzene"
    reagent = "excess phenyl magnesium bromide"
    print(f"Starting with the reaction between {reactant} and {reagent}.")
    print("The reaction proceeds via nucleophilic substitution on the aromatic ring.")
    
    # Halogen reactivity order: I > Br > Cl > F
    print("\nStep 1: The reactivity of the halogens on the benzene ring is determined by the C-X bond strength.")
    print("The order of reactivity is I > Br. The C-I bond is the weakest and will react first.")
    
    intermediate_product = "1,3-dibromo-2-phenylbenzene"
    print(f"--> The iodine at position 2 is replaced by a phenyl group, forming: {intermediate_product}")

    print("\nStep 2: The reaction conditions specify 'excess' reagent and 'reflux' (heating).")
    print("These conditions drive the reaction to completion, substituting the less reactive C-Br bonds as well.")
    
    final_product_description = "1,2,3-triphenylbenzene"
    print(f"--> The two bromine atoms at positions 1 and 3 are also replaced by phenyl groups.")
    print(f"--> The final organic product is: {final_product_description}")

    # --- IUPAC Name Construction ---
    print("\n--- Constructing the Final IUPAC Name ---")
    
    # Define the components of the IUPAC name
    locants = [1, 2, 3]
    prefix = "tri"
    substituent = "phenyl"
    parent = "benzene"
    
    # As requested, outputting each number in the name's construction
    num1, num2, num3 = locants
    print(f"The positions of the substituents (locants) are {num1}, {num2}, and {num3}.")
    
    # Constructing and printing the final name
    final_name = f"{num1},{num2},{num3}-{prefix}{substituent}{parent}"
    
    print("\nThe final IUPAC name is:")
    print(final_name)

# Execute the function
solve_reaction()