import re

def solve_wittig_reaction():
    """
    This script identifies the product of a Wittig reaction, provides the full
    reaction equation with IUPAC names, and extracts all numerical locants
    from the names in the equation as requested.
    """
    # Define the IUPAC names of the reactants and products
    aldehyde_name = "2,2-dimethylpropanal"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    alkene_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct_name = "triphenylphosphine oxide"

    # Explain the process
    print("The Wittig reaction involves an aldehyde and a ylide.")
    print(f"Reactant 1 (aldehyde): Pivalaldehyde, or {aldehyde_name}")
    print(f"Reactant 2 (ylide): {ylide_name}")
    print("\nThe carbonyl oxygen from the aldehyde and the PPh3 group from the ylide are eliminated to form a new C=C double bond.")
    print("The two products are the desired alkene and triphenylphosphine oxide.\n")

    # Construct and print the final chemical equation using IUPAC names
    final_equation = f"{aldehyde_name} + {ylide_name} -> {alkene_product_name} + {byproduct_name}"
    print("Final Reaction Equation:")
    print(final_equation)

    # Find and print all numbers from the IUPAC names in the final equation
    print("\nThe numbers from the IUPAC names in the final equation are:")
    numbers_in_equation = re.findall(r'\d+', final_equation)
    for number in numbers_in_equation:
        print(number)

solve_wittig_reaction()