# You may need to install the rdkit library first:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem

def analyze_molecules(smiles1, smiles2):
    """
    Analyzes the relationship between two molecules provided as SMILES strings.
    Determines if they are identical, constitutional isomers, or different compounds.
    """
    # Create molecule objects from SMILES strings
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if not mol1 or not mol2:
        print("Error: Invalid SMILES string provided.")
        return

    # Calculate molecular formulas
    formula1 = AllChem.CalcMolFormula(mol1)
    formula2 = AllChem.CalcMolFormula(mol2)

    print(f"Molecule 1: {smiles1}")
    print(f"Molecule 2: {smiles2}")
    print(f"Formula 1: {formula1}")
    print(f"Formula 2: {formula2}")
    print("-" * 20)

    # Compare formulas
    if formula1 != formula2:
        print("Result: The molecules have different molecular formulas. They are different compounds.")
        return

    # If formulas are the same, they could be isomers or identical.
    # To check for constitutional isomers vs. identical molecules, we can compare their canonical SMILES.
    # Canonical SMILES is a standardized representation. If they are the same, the molecules have the same connectivity.
    canonical_smiles1 = Chem.MolToSmiles(mol1, canonical=True)
    canonical_smiles2 = Chem.MolToSmiles(mol2, canonical=True)

    if canonical_smiles1 == canonical_smiles2:
        print("Result: The molecules have the same formula and the same connectivity. They are IDENTICAL.")
        print("Note: This does not distinguish between stereoisomers like enantiomers, which would require 3D analysis.")
    else:
        print("Result: The molecules have the same formula but different connectivity. They are CONSTITUTIONAL ISOMERS.")


# --- Example Usage ---
# Example 1: Constitutional Isomers (Ethanol vs. Dimethyl Ether)
print("--- Analysis 1: Ethanol vs. Dimethyl Ether ---")
ethanol_smiles = "CCO"
dimethyl_ether_smiles = "COC"
analyze_molecules(ethanol_smiles, dimethyl_ether_smiles)

print("\n" + "="*40 + "\n")

# Example 2: Identical Molecules (different representations of butane)
print("--- Analysis 2: Butane vs. Butane ---")
butane1_smiles = "CCCC"
butane2_smiles = "C(C)CC" # A different but valid way to write butane
analyze_molecules(butane1_smiles, butane2_smiles)
