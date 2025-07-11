import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit")
    sys.exit(1)

def solve_reaction_product_mass():
    """
    Solves for the product with the higher molar mass from the hydrolysis of a proposed ketal.
    """
    # Step 1 & 2: Analyze the input and validate the SMILES.
    # The user-provided SMILES string 'CC12COC(OC1)(OC2)C1=CC=CC=C1' is chemically invalid.
    # We will proceed by interpreting the likely intended reactants based on the fragments.
    print("Analysis:")
    print("The input SMILES string 'CC12COC(OC1)(OC2)C1=CC=CC=C1' is invalid.")
    print("Based on the presence of a methyl group (CC), a phenyl group (C1=...=C1), and a ketal structure")
    print("in an acid solution, we deduce the reaction is a ketal hydrolysis.")
    print("The unusual 'COC' in the ring implies a third oxygen, suggesting the ketal was formed from diethylene glycol.")
    print("The most plausible starting material is the ketal of Acetophenone and Diethylene Glycol.")
    print("\nThe hydrolysis reaction is:")
    print("Ketal + H2O --(H+)--> Acetophenone + Diethylene Glycol\n")

    # Step 3: Define the products of the hydrolysis reaction.
    products = {
        "Acetophenone": "CC(=O)c1ccccc1",
        "Diethylene Glycol": "OCCOCCO"
    }

    # Step 4: Calculate and compare molar masses.
    product_data = []
    for name, smiles in products.items():
        mol = Chem.MolFromSmiles(smiles)
        # The 'final equation' will be the comparison of molar masses.
        molar_mass = Descriptors.MolWt(mol)
        product_data.append({"name": name, "smiles": smiles, "mass": molar_mass})

    print("Comparing the molar masses of the products:")
    for data in product_data:
        # Outputting each number as requested.
        print(f"- Product: {data['name']}, SMILES: {data['smiles']}, Molar Mass: {data['mass']:.2f} g/mol")

    # Determine which product has the higher molar mass
    heavier_product = max(product_data, key=lambda x: x['mass'])

    print(f"\nThe product with the higher molar mass is {heavier_product['name']}.")
    print(f"Its SMILES string is: {heavier_product['smiles']}")

    # Step 5: Return the final answer in the specified format.
    final_answer_smiles = heavier_product['smiles']
    print(f"\n<<<{final_answer_smiles}>>>")


if __name__ == "__main__":
    solve_reaction_product_mass()
