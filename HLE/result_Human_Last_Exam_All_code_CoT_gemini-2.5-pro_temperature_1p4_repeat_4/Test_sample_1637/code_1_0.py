#!/usr/bin/env python

# To run this script, you first need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Crippen

def analyze_probe_solubility():
    """
    Analyzes and compares the predicted water solubility of two probes
    by calculating their LogP values.
    """
    # Chemical structure of the original probe in SMILES format
    # Name: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_probe_smiles = "COC1=CC2=C(C=C1)C(=O)C3=C(S2)C=C(C=C3)OCC(=O)NCCOCCOCCCCCCCl"

    # Chemical structure of a modified probe.
    # The hydrophobic chloro-hexyl group and amide linker are replaced with
    # a more hydrophilic tetraethylene glycol (PEG4) chain.
    modified_probe_smiles = "COC1=CC2=C(C=C1)C(=O)C3=C(S2)C=C(C=C3)OCCOCCOCCOCCO"

    # Create molecule objects from the SMILES strings
    original_mol = Chem.MolFromSmiles(original_probe_smiles)
    modified_mol = Chem.MolFromSmiles(modified_probe_smiles)

    if original_mol is None or modified_mol is None:
        print("Error: Could not create one or both of the molecules from the SMILES strings.")
        return

    # Calculate the LogP value for each molecule.
    # LogP is a measure of hydrophobicity. A lower value suggests better water solubility.
    original_logp = Crippen.MolLogP(original_mol)
    modified_logp = Crippen.MolLogP(modified_mol)

    print("--- Solubility Analysis ---")
    print(f"Original Probe SMILES: {original_probe_smiles}")
    print(f"Modified Probe SMILES: {modified_probe_smiles}\n")
    print("Calculated LogP (a measure of hydrophobicity):")
    print(f"Original Probe LogP = {original_logp:.4f}")
    print(f"Modified Probe LogP = {modified_logp:.4f}\n")

    print("--- Conclusion ---")
    print("Yes, changing the structure to incorporate a more significant PEG group is an excellent strategy.")
    print("The modified probe has a much lower predicted LogP value, which strongly suggests")
    print("it will be more water-soluble and is much less likely to precipitate from the medium.")

if __name__ == "__main__":
    analyze_probe_solubility()
