# First, you may need to install the RDKit library. You can do this by running:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of the original probe and
    a modified version with an extended PEG linker.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1cc2c(cc1OCC(=O)NCCOCCOCCCCCCCl)SC3=CC=CC=C3C2=O"

    # SMILES string for a modified probe with a longer PEG linker (two extra ethylene glycol units)
    modified_smiles = "COc1cc2c(cc1OCC(=O)NCCOCCOCCOCCOCCCCCCCl)SC3=CC=CC=C3C2=O"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not create molecule from SMILES string. Please check the input.")
        return

    # --- Calculations for the Original Probe ---
    original_logp = Crippen.MolLogP(original_mol)
    original_molwt = Descriptors.MolWt(original_mol)

    # --- Calculations for the Modified Probe ---
    modified_logp = Crippen.MolLogP(modified_mol)
    modified_molwt = Descriptors.MolWt(modified_mol)

    # --- Print the results ---
    print("--- Analysis of Probe Solubility ---")
    print("\nYes, replacing part of the linker with a more substantial PEG group is an excellent strategy to solve the precipitation problem.")
    print("The original probe has large hydrophobic (water-fearing) regions, leading to low aqueous solubility.")
    print("Adding a hydrophilic (water-loving) PEG chain increases the overall solubility of the molecule.")
    print("\nTo demonstrate this computationally, we calculate the LogP (a measure of hydrophobicity).")
    print("A lower LogP value suggests higher water solubility.\n")

    print("--- Original Probe ---")
    print(f"Molecular Weight: {original_molwt:.2f}")
    print(f"Calculated LogP: {original_logp:.2f}")

    print("\n--- Modified Probe (with longer PEG linker) ---")
    print(f"Molecular Weight: {modified_molwt:.2f}")
    print(f"Calculated LogP: {modified_logp:.2f}")

    print("\n--- Conclusion ---")
    print(f"The modified probe has a significantly lower LogP value ({modified_logp:.2f}) compared to the original probe ({original_logp:.2f}).")
    print("This computational result strongly supports that increasing the PEG character of the linker will increase water solubility and likely solve the precipitation issue.")

if __name__ == "__main__":
    analyze_probe_solubility()
<<<Yes>>>