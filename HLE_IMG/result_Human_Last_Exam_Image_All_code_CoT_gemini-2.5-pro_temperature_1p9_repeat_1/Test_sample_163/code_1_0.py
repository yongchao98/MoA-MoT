# To run this code, you may need to install the rdkit library:
# pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def identify_reaction_products():
    """
    Identifies and characterizes the two major products (A and B) of the given reaction.
    """
    # The reaction of styrene with tert-butyl peroxybenzoate adds a benzoyloxy group (-OCOPh)
    # and a tert-butoxy group (-OtBu) across the double bond. The two major products are regioisomers.

    # Product A: benzoyloxy group on the benzylic carbon, tert-butoxy on the terminal carbon.
    # Name: 2-(tert-butoxy)-1-phenylethyl benzoate
    smiles_A = "c1ccc(C(OC(=O)c2ccccc2)COC(C)(C)C)cc1"

    # Product B: tert-butoxy group on the benzylic carbon, benzoyloxy on the terminal carbon.
    # Name: 2-(tert-butoxy)-2-phenylethyl benzoate
    smiles_B = "c1ccc(C(OC(C)(C)C)COC(=O)c2ccccc2)cc1"

    mol_A = Chem.MolFromSmiles(smiles_A)
    mol_B = Chem.MolFromSmiles(smiles_B)

    print("The reaction produces two major products, A and B, which are regioisomers.\n")

    # --- Product A ---
    print("Product A:")
    print("Name: 2-(tert-butoxy)-1-phenylethyl benzoate")
    print("Structure: The benzoyloxy group is attached to the carbon bearing the phenyl ring.")
    print(f"SMILES: {smiles_A}")
    print(f"Molecular Formula: {CalcMolFormula(mol_A)}")
    print(f"Molecular Weight: {Descriptors.MolWt(mol_A):.4f}")
    print("-" * 20)

    # --- Product B ---
    print("Product B:")
    print("Name: 2-(tert-butoxy)-2-phenylethyl benzoate")
    print("Structure: The tert-butoxy group is attached to the carbon bearing the phenyl ring.")
    print(f"SMILES: {smiles_B}")
    print(f"Molecular Formula: {CalcMolFormula(mol_B)}")
    print(f"Molecular Weight: {Descriptors.MolWt(mol_B):.4f}")
    print("-" * 20)


identify_reaction_products()
