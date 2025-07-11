import re

def solve_wittig_reaction():
    """
    This function analyzes a Wittig reaction and prints the name of the product.
    """
    # 1. Define reactants
    aldehyde = "pivalaldehyde"
    ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    byproduct = "triphenylphosphine oxide"

    # 2. Determine the product structure and name
    # Pivalaldehyde: (CH3)3C-CHO
    # Ylide: Ph3P=CH-CH2-(C6H4Cl at position 2)
    # The reaction replaces C=O with C=C from the ylide.
    # Product structure: (CH3)3C-CH=CH-CH2-(2-chlorophenyl)
    product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 3. Print the reaction summary and the result
    print("The Wittig reaction converts an aldehyde/ketone and a phosphonium ylide into an alkene.")
    print("-" * 70)
    print(f"Aldehyde Reactant: {aldehyde}")
    print(f"Ylide Reactant: {ylide}")
    print("\nReaction Equation:")
    print(f"{aldehyde} + {ylide} -> {product_name} + {byproduct}")
    print("-" * 70)
    print(f"The IUPAC name of the major organic product is:")
    print(product_name)
    print("\nAs requested, here are the numbers present in the product's IUPAC name:")
    
    # 4. Extract and print each number from the name
    numbers = re.findall(r'\d', product_name)
    for num in numbers:
        print(num)

if __name__ == "__main__":
    solve_wittig_reaction()