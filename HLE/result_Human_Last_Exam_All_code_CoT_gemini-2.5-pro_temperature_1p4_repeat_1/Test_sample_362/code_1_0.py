import re

def wittig_reaction_product():
    """
    This function identifies the products of the specified Wittig reaction,
    prints the reaction equation, and extracts all numbers from the chemical formulas.
    """
    # Reactant names and formulas
    aldehyde_name = "Pivalaldehyde"
    aldehyde_formula = "C5H10O"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_formula = "C26H22ClP"

    # Product names and formulas
    alkene_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    alkene_product_formula = "C13H17Cl"
    byproduct_name = "Triphenylphosphine oxide"
    byproduct_formula = "C18H15PO"

    # Print the product names
    print(f"The Wittig reaction produces two compounds:")
    print(f"1. Main Product: {alkene_product_name} ({alkene_product_formula})")
    print(f"2. Byproduct: {byproduct_name} ({byproduct_formula})")
    print("-" * 20)

    # Construct and print the full reaction equation
    reaction_equation = (
        f"{aldehyde_formula} ({aldehyde_name}) + "
        f"{ylide_formula} (Ylide) -> "
        f"{alkene_product_formula} (Alkene) + "
        f"{byproduct_formula} (Byproduct)"
    )
    print("Full Reaction Equation:")
    print(reaction_equation)
    print("-" * 20)

    # As requested, extract and print all numbers from the chemical formulas in the equation
    # The equation for formula extraction is: C5H10O + C26H22ClP -> C13H17Cl + C18H15PO
    equation_formulas_only = f"{aldehyde_formula} + {ylide_formula} -> {alkene_product_formula} + {byproduct_formula}"
    numbers = re.findall(r'\d+', equation_formulas_only)
    
    print("The numbers in the final equation's formulas are:")
    # Print each number as requested
    for num in numbers:
        print(num)

# Execute the function
wittig_reaction_product()