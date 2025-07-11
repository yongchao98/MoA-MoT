import sys

def chemical_reaction_analysis():
    """
    Analyzes the Diels-Alder reaction between 1,3-butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the result.
    """
    print("--- Reaction Analysis ---")
    print("Reactant 1: 1,3-Butadiene (a conjugated diene)")
    print("Reactant 2: 1,1-dichloro-2,2-difluoroethene (a dienophile)")
    print("Reaction Type: Diels-Alder [4+2] Cycloaddition")
    print("\n--- Mechanism ---")
    print("The four pi-electrons from the diene and two pi-electrons from the dienophile")
    print("rearrange in a concerted step to form a new six-membered ring (a cyclohexene derivative).")
    print("\n--- Product Formation ---")
    print("The atoms from 1,3-butadiene form four carbons of the new ring, including the new double bond.")
    print("The two carbons from 1,1-dichloro-2,2-difluoroethene form the other two carbons of the ring.")
    print("The chlorine and fluorine atoms remain attached to their respective carbons.")
    
    # Define reactants and product
    reactant1_name = "1,3-butadiene"
    reactant2_name = "1,1-dichloro-2,2-difluoroethene"
    product_name = "3,3-dichloro-4,4-difluorocyclohexene"

    print("\n--- Final Equation ---")
    # As requested, outputting the names which include the relevant numbers
    print(f"The reaction is: {reactant1_name} + {reactant2_name} -> {product_name}")
    print("\nThe product of the reaction is 3,3-dichloro-4,4-difluorocyclohexene.")

if __name__ == '__main__':
    chemical_reaction_analysis()