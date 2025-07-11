import re

def solve_chemistry_problem():
    """
    This function identifies the final product of the described chemical synthesis
    and prints its name and chemical formula.
    """
    # Based on the step-by-step analysis of the reaction sequence,
    # the final product, Compound 4, is identified.
    final_product_name = "Phthalide"
    chemical_formula = "C8H6O2"

    print(f"The final product, Compound 4, is: {final_product_name}")
    print(f"The chemical formula is: {chemical_formula}")

    # The prompt asks to output each number in the final equation (formula).
    numbers = re.findall(r'\d+', chemical_formula)
    print("The numbers in the final chemical formula are:")
    for num in numbers:
        print(num)

solve_chemistry_problem()