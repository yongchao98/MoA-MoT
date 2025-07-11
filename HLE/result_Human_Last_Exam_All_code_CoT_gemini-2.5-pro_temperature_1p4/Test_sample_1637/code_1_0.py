# To run this code, you first need to install the rdkit library.
# You can do this by running the following command in your terminal or command prompt:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of the original probe and a
    modified version with a longer PEG linker.
    """
    # SMILES is a string representation of a chemical structure.
    # Original probe: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    # This molecule has a linker with two ethylene glycol units.
    original_smiles = "COc1ccc2c(c1)SC1=C(C2=O)C=C(OCC(=O)NCCOCCOCCCCCCCl)C=C1"

    # Modified probe with a longer PEG linker.
    # We extend the linker from 2 ethylene glycol units to 4.
    # This corresponds to: N-(2-(2-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethoxy)ethoxy)ethyl)-...
    modified_smiles = "COc1ccc2c(c1)SC1=C(C2=O)C=C(OCC(=O)NCCOCCOCCOCCOCCCCCCCl)C=C1"

    # Create molecule objects from the SMILES strings
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    # Calculate cLogP (a measure of hydrophobicity) for both molecules
    # A lower cLogP value predicts better water solubility.
    original_clogp = Descriptors.MolLogP(original_mol)
    modified_clogp = Descriptors.MolLogP(modified_mol)

    print("Analysis of Probe Solubility:")
    print("-" * 30)
    print(f"Original Probe SMILES: {original_smiles}")
    print(f"Modified Probe SMILES: {modified_smiles}")
    print("-" * 30)
    
    # Print the final equation comparing the two values
    print("Predicted Hydrophobicity (cLogP):")
    print(f"Original Probe cLogP = {original_clogp:.2f}")
    print(f"Modified Probe cLogP = {modified_clogp:.2f}")
    print("-" * 30)
    
    if modified_clogp < original_clogp:
        print("Conclusion: The modified probe with the longer PEG linker has a lower cLogP.")
        print("This predicts a significant improvement in aqueous solubility, which should help solve the precipitation problem.")
    else:
        print("Conclusion: The modification did not lower the cLogP as expected.")

if __name__ == '__main__':
    analyze_probe_solubility()
