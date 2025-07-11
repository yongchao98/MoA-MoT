import re

def solve_reaction():
    """
    Analyzes the Diels-Alder reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the result.
    """
    print("Step 1: Identifying the reactants and reaction type.")
    print("Reactant 1 (Diene): 1,3-Butadiene (C4H6)")
    print("Reactant 2 (Dienophile): 1,1-dichloro-2,2-difluoroethene (C2Cl2F2)")
    print("This is a [4+2] cycloaddition known as the Diels-Alder reaction, which forms a six-membered ring.\n")

    print("Step 2: Predicting the product structure.")
    print("The 4 carbons of the diene and the 2 carbons of the dienophile combine.")
    print("A new cyclohexene ring is formed. The substituents (Cl and F) from the dienophile remain attached to their carbons.\n")

    print("Step 3: Naming the product.")
    product_name = "3,3-dichloro-4,4-difluorocyclohex-1-ene"
    print(f"The IUPAC name of the product is: {product_name}\n")

    print("Step 4: The final chemical equation.")
    equation = "C4H6 + C2Cl2F2 -> C6H6Cl2F2"
    print(f"The balanced chemical equation is: {equation}\n")

    print("Step 5: Extracting numbers from the final equation as requested.")
    # Use regex to find all numbers in the equation string
    numbers = re.findall(r'\d+', equation)
    # The prompt asks to output each number in the final equation.
    # For C4H6 + C2Cl2F2 -> C6H6Cl2F2, the numbers are 4, 6, 2, 2, 2, 6, 6, 2, 2.
    # Note: C2Cl2F2 implies C2Cl2F2, so the numbers are 2, 2, 2.
    # The product C6H6Cl2F2 implies C6H6Cl2F2, so the numbers are 6, 6, 2, 2.
    print("The numbers present in the chemical formulas of the equation are:")
    print(', '.join(numbers))


solve_reaction()
<<<3,3-dichloro-4,4-difluorocyclohex-1-ene>>>