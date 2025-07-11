import re

def solve_wittig_reaction():
    """
    This function analyzes a Wittig reaction and prints the results.
    """
    # Step 1: Define the reactants
    aldehyde_name = "Pivalaldehyde (2,2-dimethylpropanal)"
    aldehyde_structure = "(CH3)3C-CHO"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_structure = "(2-Cl-C6H4)-CH2-CH=P(Ph)3"

    # Step 2: Describe the reaction
    print("### Wittig Reaction Analysis ###\n")
    print("The Wittig reaction creates an alkene from a carbonyl compound and a phosphorus ylide.\n")
    print("Reactants:")
    print(f"1. Aldehyde: {aldehyde_name}")
    print(f"   Structure: {aldehyde_structure}")
    print(f"2. Ylide (Wittig Reagent): {ylide_name}")
    print(f"   Structure: {ylide_structure}\n")

    # Step 3: Explain the transformation
    print("Mechanism:")
    print("The oxygen atom of the aldehyde's carbonyl group (C=O) is swapped with the")
    print("carbon group from the ylide (C=P) to form a new carbon-carbon double bond (C=C).\n")

    # Step 4: Identify the products
    product_name_long = "(E)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    product_structure = "(2-Cl-C6H4)-CH2-CH=CH-C(CH3)3"
    byproduct_name = "Triphenylphosphine oxide"
    byproduct_structure = "(Ph)3P=O"
    
    print("Products:")
    print(f"1. Main Product: {product_name_long}")
    print(f"   Structure: {product_structure}")
    print("   Stereochemistry: The reaction involves a bulky aldehyde (pivalaldehyde) and a")
    print("   non-stabilized ylide. The steric hindrance favors the formation of the more")
    print("   stable E (trans) isomer as the major product.\n")
    print(f"2. Byproduct: {byproduct_name}")
    print(f"   Structure: {byproduct_structure}\n")

    # Step 5: Fulfill the request to show the final equation and numbers
    print("---")
    print("Final Reaction Equation:")
    print(f"{aldehyde_name} + {ylide_name} ---> {product_name_long} + {byproduct_name}\n")
    
    # Extract numbers from the IUPAC name
    numbers_in_name = re.findall(r'\d+', product_name_long)
    
    print("As requested, here are the numbers present in the IUPAC name of the main product in the final equation:")
    print(f"Product Name: '{product_name_long}'")
    print("Numbers:", ', '.join(numbers_in_name))
    print("---")

solve_wittig_reaction()
