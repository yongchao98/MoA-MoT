import re

def solve_reaction_products():
    """
    Identifies the products of the reaction between styrene and tert-butyl peroxybenzoate,
    and prints details of the overall chemical equation.
    """
    # Define product structures based on chemical principles.
    # The reaction adds a benzoyloxy group (-OBz) and a tert-butoxy group (-OtBu)
    # across the double bond of styrene (Ph-CH=CH2).
    product_A_name = "1-(Benzoyloxy)-2-(tert-butoxy)-1-phenylethane"
    product_A_structure = "C6H5-CH(OCOC6H5)-CH2-OC(CH3)3"

    product_B_name = "2-(Benzoyloxy)-1-(tert-butoxy)-1-phenylethane"
    product_B_structure = "C6H5-CH(OC(CH3)3)-CH2-OCOC6H5"

    print("The two major products, A and B, are the following constitutional isomers:\n")
    print(f"Product A: {product_A_name}")
    print(f"           Structure: {product_A_structure}\n")
    print(f"Product B: {product_B_name}")
    print(f"           Structure: {product_B_structure}\n")

    # Define molecular formulas and stoichiometry for the overall equation
    reactants = {'Styrene': 'C8H8', 'tert-Butyl peroxybenzoate': 'C11H14O3'}
    products_formula = 'C19H22O3'
    stoichiometry = 1

    # Construct and print the final equation
    equation = (
        f"{stoichiometry} {reactants['Styrene']} + "
        f"{stoichiometry} {reactants['tert-Butyl peroxybenzoate']} -> "
        f"{stoichiometry} {products_formula}"
    )
    print("The overall balanced chemical equation is:")
    print(equation)
    print("\nAs requested, here are all the numbers from the final equation:")

    # Extract and print numbers from the equation string
    all_numbers = re.findall(r'\d+', equation)
    for number in all_numbers:
        print(number)

# Run the analysis
solve_reaction_products()