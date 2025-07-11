import re

def identify_compound_3():
    """
    This function outlines the synthesis and identifies the final product, Compound 3.
    The reaction sequence is as follows:
    1. Terpinolene is epoxidized at its exocyclic double bond with m-CPBA to give compound 1 (a spiro-epoxide).
    2. Compound 1 is treated with N,N-dimethyl thioformamide and acid to convert the epoxide into a thiirane (compound 2).
    3. Compound 2 is reduced with LiAlH4. The hydride attacks the less sterically hindered spiro-carbon of the thiirane,
       leading to a reductive ring-opening to form the final product, compound 3 (a thiol).
    """

    # Identifiers for the final product, Compound 3
    iupac_name = "4-(2-sulfanylpropan-2-yl)-1-methylcyclohex-1-ene"
    molecular_formula = "C10H18S"
    smiles_string = "CC1=CCC(C(C)(C)S)CC1"

    # Per the instructions, extract and print the numbers from the IUPAC name
    # This fulfills the "output each number in the final equation!" requirement.
    numbers_in_name = re.findall(r'\d+', iupac_name)

    print("Identification of Final Product (Compound 3):")
    print("-" * 40)
    print(f"IUPAC Name: {iupac_name}")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"SMILES String: {smiles_string}")
    print("-" * 40)
    print("Numbers from the IUPAC Name:")
    for number in numbers_in_name:
        print(number)

if __name__ == "__main__":
    identify_compound_3()
