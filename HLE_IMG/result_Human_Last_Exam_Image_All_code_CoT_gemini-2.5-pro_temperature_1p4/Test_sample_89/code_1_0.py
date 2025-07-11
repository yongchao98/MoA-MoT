import sys
try:
    from rdkit import Chem
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit()

def get_product_structures():
    """
    This function provides the SMILES strings for the three products A, B, and C,
    based on the analysis of the reaction description.
    """
    # Based on the text, Product A is described as a "primary amide".
    # A plausible structure derived from the proline part of the starting material
    # is N-acetylpyrrolidine-2-carboxamide.
    smiles_A = "CC(=O)N1CCCC1C(=O)N"

    # Product B is described as a "tethered imide" formed from a Michael addition.
    # It contains the dihydropyrrole "tether", an acetyl group, and a group
    # derived from the methyl propiolate adduct.
    # SMILES for N-acetyl-N-(3-(methoxycarbonyl)acryloyl)-(4,5-dihydro-3H-pyrrol-2-yl)amine
    smiles_B = "COC(=O)/C=C/C(=O)N(C1=NCCC1)C(=O)C"

    # Product C is described as a "tethered imide" from reaction with acetic anhydride.
    # The most logical structure is the di-acetylated dihydropyrrole amine.
    smiles_C = "CC(=O)N(C1=NCCC1)C(=O)C"

    return smiles_A, smiles_B, smiles_C

def main():
    """
    Main function to print the structures of the products.
    """
    smiles_A, smiles_B, smiles_C = get_product_structures()

    print("The proposed structures for products A, B, and C are provided below in SMILES format.")
    print("-" * 40)
    print("Structure of Product A:")
    print(smiles_A)
    # Validate SMILES and print formula as a check
    mol_A = Chem.MolFromSmiles(smiles_A)
    formula_A = Chem.rdMolDescriptors.CalcMolFormula(mol_A)
    print(f"Molecular Formula: {formula_A}")
    print("-" * 40)

    print("Structure of Product B:")
    print(smiles_B)
    mol_B = Chem.MolFromSmiles(smiles_B)
    formula_B = Chem.rdMolDescriptors.CalcMolFormula(mol_B)
    print(f"Molecular Formula: {formula_B}")
    print("-" * 40)

    print("Structure of Product C:")
    print(smiles_C)
    mol_C = Chem.MolFromSmiles(smiles_C)
    formula_C = Chem.rdMolDescriptors.CalcMolFormula(mol_C)
    print(f"Molecular Formula: {formula_C}")
    print("-" * 40)


if __name__ == "__main__":
    main()