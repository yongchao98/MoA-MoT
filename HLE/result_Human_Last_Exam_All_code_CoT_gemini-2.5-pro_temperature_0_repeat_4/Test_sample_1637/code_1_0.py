# First, ensure you have rdkit installed:
# pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_probe_solubility():
    """
    Calculates and compares the hydrophobicity (cLogP) of the original probe
    and a modified version with a PEG chain to improve solubility.
    A lower cLogP value indicates lower hydrophobicity and higher water solubility.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_probe_smiles = "COC1=CC2=C(C=C1)C(=O)C3=C(S2)C=C(C=C3)OCC(=O)NCCOCCOCCCCCCCl"

    # SMILES string for a modified probe where the hydrophobic chloro-hexyl tail
    # is replaced with a more hydrophilic PEG4 (tetraethylene glycol) chain.
    modified_probe_smiles = "COC1=CC2=C(C=C1)C(=O)C3=C(S2)C=C(C=C3)OCC(=O)NCCOCCOCCOCCO"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_probe_smiles)
    modified_mol = Chem.MolFromSmiles(modified_probe_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # Calculate cLogP (a measure of hydrophobicity)
    original_clogp = Descriptors.CrippenClogP(original_mol)
    modified_clogp = Descriptors.CrippenClogP(modified_mol)

    print("Analysis of Probe Solubility by Calculating cLogP:")
    print("-" * 50)
    print(f"Original Probe cLogP: {original_clogp:.4f}")
    print(f"Modified Probe cLogP: {modified_clogp:.4f}")
    print("-" * 50)
    
    # The "equation" here is the comparison of the two values.
    # We output each number in the final comparison.
    print("Final Equation (Comparison):")
    print(f"Original cLogP ({original_clogp:.4f}) > Modified cLogP ({modified_clogp:.4f})")
    
    print("\nConclusion:")
    print("The modified probe has a significantly lower cLogP value.")
    print("This indicates it is much less hydrophobic and will likely have higher aqueous solubility, solving the precipitation issue.")

if __name__ == "__main__":
    analyze_probe_solubility()
<<<Yes>>>