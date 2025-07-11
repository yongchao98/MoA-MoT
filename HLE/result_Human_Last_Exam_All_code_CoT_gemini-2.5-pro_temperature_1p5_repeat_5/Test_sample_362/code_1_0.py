import re

def solve_wittig_reaction():
    """
    This function identifies the product of a specific Wittig reaction,
    prints the reaction equation, and highlights the numbering in the product's IUPAC name.
    """
    
    # 1. Identify the structure and name of the final product.
    # Reactants: Pivalaldehyde ((CH3)3C-CHO) and the ylide (Ph3P=CH-CH2-(o-Cl-C6H4))
    # The reaction swaps the =O of the aldehyde with the =C< part of the ylide.
    # Resulting alkene: (CH3)3C-CH=CH-CH2-(o-Cl-C6H4)
    # The non-stabilized ylide favors the Z-isomer.
    product_iupac_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    product_structure = "(CH3)3C-CH=CH-CH2-(2-Cl-C6H4)"

    # 2. Define the overall reaction equation.
    # The task requires outputting the numbers "in the final equation".
    # This refers to the locant numbers that define the structure of the product.
    final_equation = f"Pivalaldehyde + Ylide  ->  {product_structure} + Triphenylphosphine oxide"

    # 3. Extract the numbers from the IUPAC name.
    numbers_in_name = re.findall(r'\d', product_iupac_name)

    # 4. Print the results clearly.
    print("The primary product of the Wittig reaction is an alkene.")
    print("\n--- Final Reaction Equation ---")
    print(final_equation)
    
    print("\nThe IUPAC name of the major organic product is:")
    print(product_iupac_name)

    print("\nThe numbers from the product's name, which define the final structure, are:")
    # Print each number separated by a space
    print(' '.join(numbers_in_name))

solve_wittig_reaction()