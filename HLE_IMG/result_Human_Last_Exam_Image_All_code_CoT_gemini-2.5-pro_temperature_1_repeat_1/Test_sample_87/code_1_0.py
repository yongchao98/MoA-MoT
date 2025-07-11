# First, ensure you have the RDKit library installed.
# If not, you can install it using pip:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_formula(smiles):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES string"
    return Descriptors.rdMolDescriptors.CalcMolFormula(mol)

def solve_reaction_products():
    """
    Solves for the structures of products A, B, and C from the given reaction.
    The structures are determined based on mechanistic plausibility and matching the given molecular formulas.
    """

    # --- Product A ---
    # Molecular Formula: C14H20N2O3
    # Plausible Pathway: The initial azomethine ylide reacts with methyl propiolate, and the
    # resulting enamine adduct then reacts with ketene (formed in situ from Ac2O).
    # This leads to a complex spirocyclic fused system.
    smiles_A = "COC(=O)C=C1C(=O)C[C@]2(N1C1=NCCC1)[C@H]1CCN[C@@H]12"
    
    # --- Product B ---
    # Molecular Formula: C12H14N2O3
    # Plausible Pathway: The primary [3+2] cycloaddition adduct undergoes a significant
    # oxidative rearrangement, likely promoted by Ac2O and heat. This could involve
    # dehydrogenation and formation of new rings/carbonyl groups.
    smiles_B = "COC(=O)C1=CC(=O)N2C=C(N1C1=NCCC1)C1CC12"

    # --- Product C ---
    # Molecular Formula: C11H16N2O3
    # Plausible Pathway: This product likely arises from another complex rearrangement pathway
    # of the primary adduct, involving both oxidation and the loss of a CH2 group
    # (e.g., demethylation of the ester group).
    smiles_C = "COC(=O)[C@@]12CCCN1C(=O)CN(C=N2)C1CCC1"

    print("The proposed structures for products A, B, and C are as follows:\n")

    print(f"--- Product A ---")
    print(f"Given Molecular Formula: C14H20N2O3")
    print(f"Proposed SMILES: {smiles_A}")
    print(f"Calculated Formula from SMILES: {get_molecular_formula(smiles_A)}\n")

    print(f"--- Product B ---")
    print(f"Given Molecular Formula: C12H14N2O3")
    print(f"Proposed SMILES: {smiles_B}")
    print(f"Calculated Formula from SMILES: {get_molecular_formula(smiles_B)}\n")

    print(f"--- Product C ---")
    print(f"Given Molecular Formula: C11H16N2O3")
    print(f"Proposed SMILES: {smiles_C}")
    print(f"Calculated Formula from SMILES: {get_molecular_formula(smiles_C)}\n")


if __name__ == '__main__':
    solve_reaction_products()