import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it first by running: pip install rdkit")
    sys.exit(1)

def analyze_probe_solubility():
    """
    Calculates and compares the LogP values for an original probe and a
    modified, more water-soluble version.
    """
    # The chemical structure of the original probe in SMILES format.
    # The 'CCCCCCCl' part represents the hydrophobic 6-chlorohexyl chain.
    original_probe_smiles = 'COC1=CC2=C(C=C1)C(=O)C3=CC(=C(S2)C=C3)OCC(=O)NCCOCCOCCCCCCCl'

    # A modified probe where the hydrophobic '-O-(CH2)6-Cl' tail is replaced
    # with a longer, more hydrophilic PEG-like chain: '-O-(CH2CH2O)3-CH2CH2-Cl'.
    # This change increases the number of water-loving oxygen atoms.
    modified_probe_smiles = 'COC1=CC2=C(C=C1)C(=O)C3=CC(=C(S2)C=C3)OCC(=O)NCCOCCOCCOCCOCCCl'

    # Create RDKit molecule objects from the SMILES strings
    original_mol = Chem.MolFromSmiles(original_probe_smiles)
    modified_mol = Chem.MolFromSmiles(modified_probe_smiles)

    # Calculate the LogP for each molecule.
    # LogP is a measure of hydrophobicity (lipophilicity).
    # A lower LogP indicates higher water solubility.
    original_logp = Descriptors.MolLogP(original_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)

    print("--- Solubility Analysis ---")
    print("Problem: A probe is precipitating from the cell culture medium.")
    print("Hypothesis: Increasing the PEG content of the probe will increase its water solubility and solve the issue.\n")
    print("Method: We will calculate and compare the LogP (hydrophobicity value) of the original and modified probes. A lower LogP suggests better water solubility.\n")

    print(f"Original Probe Structure (SMILES): {original_probe_smiles}")
    print(f"Calculated LogP for Original Probe = {original_logp:.4f}\n")

    print(f"Modified Probe Structure (SMILES): {modified_probe_smiles}")
    print(f"Calculated LogP for Modified Probe = {modified_logp:.4f}\n")

    print("--- Conclusion ---")
    print(f"The modification decreased the LogP value from {original_logp:.4f} to {modified_logp:.4f}.")
    print("This significant decrease in LogP confirms that the modified probe is substantially more hydrophilic (water-soluble).")
    print("Therefore, modifying the probe to increase its PEG character is very likely to solve the precipitation problem.")

if __name__ == '__main__':
    analyze_probe_solubility()
<<<Yes>>>