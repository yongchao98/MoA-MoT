import re

def describe_reaction():
    """
    This script describes the product of the reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene.
    """

    # 1. Define reactants and product
    reactant1_name = "Butadiene"
    reactant1_formula = "C4H6"

    reactant2_name = "1,1-dichloro-2,2-difluoroethene"
    reactant2_formula = "C2Cl2F2"

    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    product_formula = "C6H6Cl2F2"

    # 2. Explain the reaction
    print("### Reaction Analysis ###")
    print(f"The reaction between {reactant1_name} and {reactant2_name} is a Diels-Alder reaction.")
    print("This is a [4+2] cycloaddition where the conjugated diene (butadiene) reacts with the")
    print("dienophile (1,1-dichloro-2,2-difluoroethene) to form a six-membered ring.\n")

    # 3. Show the chemical equation
    equation = f"{reactant1_formula} + {reactant2_formula} -> {product_formula}"
    print("Chemical Equation:")
    print(equation)
    print("\n")

    # 4. State the final product name
    print("Product Name:")
    print(f"The major product of this reaction is: {product_name}")
    print("\n")

    # 5. Extract and print all numbers from the equation as requested
    numbers = re.findall(r'\d', equation)
    print("Numbers from the final equation:")
    print(' '.join(numbers))


if __name__ == '__main__':
    describe_reaction()