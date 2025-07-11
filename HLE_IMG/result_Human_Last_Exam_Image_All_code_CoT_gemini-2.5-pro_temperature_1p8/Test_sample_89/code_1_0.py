from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def solve_reaction_products():
    """
    This script determines the structures of products A, B, and C based on the provided
    reaction pathways and prints their key properties. The solution is based on interpreting
    the (sometimes contradictory) text by applying known principles of organic chemistry.
    """
    # SMILES strings for the proposed structures based on the analysis
    structures = {
        "Product A: Acetamide": "CC(=O)N",
        "Product B: Hexahydropyrrolizin-3-one": "O=C1CN2C(CCC2)C1",
        "Product C: (S)-3-acetyl-pyrrolizinone derivative": "CC(=O)C1=C[C@@H]2N(C1=O)CCCC2"
    }

    print("Based on the reaction description, the proposed structures for products A, B, and C are:")
    print("-" * 70)

    for name, smiles in structures.items():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                # This check is for robustness
                print(f"Could not parse SMILES for {name}: {smiles}")
                continue
            
            # Calculate molecular formula and weight
            formula = CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)

            print(f"{name}")
            print(f"  Molecular Formula: {formula}")
            print(f"  Molecular Weight: {mw:.2f}")
            print(f"  SMILES: {smiles}")
            print("-" * 70)
        except Exception as e:
            print(f"An error occurred while processing {name}: {e}")

solve_reaction_products()