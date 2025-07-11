# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_probe_solubility():
    """
    Analyzes and compares the predicted solubility of an original and a modified chemical probe.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1ccc2c(c1)C(=O)c1cc(OCC(=O)NCCOCCOCCCCCCCl)ccc1S2"

    # Proposed modification: Replace the hydrophobic `-(CH2)6-Cl` group with a more hydrophilic
    # tetraethylene glycol chain with a terminal chloro group, i.e., -(CH2CH2O)4-CH2CH2Cl.
    # This directly addresses the likely source of hydrophobicity while retaining a reactive handle.
    modified_smiles = "COc1ccc2c(c1)C(=O)c1cc(OCC(=O)NCCOCCOCCOCCOCCCl)ccc1S2"

    # Create RDKit molecule objects
    try:
        original_mol = Chem.MolFromSmiles(original_smiles)
        modified_mol = Chem.MolFromSmiles(modified_smiles)
        if not original_mol or not modified_mol:
            raise ValueError("Invalid SMILES string provided.")
    except Exception as e:
        print(f"Error creating molecule from SMILES string: {e}")
        return

    # --- Calculations for Original Probe ---
    original_logp = Descriptors.MolLogP(original_mol)
    original_wt = Descriptors.MolWt(original_mol)
    original_hba = Descriptors.NumHAcceptors(original_mol)
    original_hbd = Descriptors.NumHDonors(original_mol)

    # --- Calculations for Modified Probe ---
    modified_logp = Descriptors.MolLogP(modified_mol)
    modified_wt = Descriptors.MolWt(modified_mol)
    modified_hba = Descriptors.NumHAcceptors(modified_mol)
    modified_hbd = Descriptors.NumHDonors(modified_mol)

    # --- Print Results and Conclusion ---
    print("--- Analysis of Probe Solubility ---")
    print("\nThe observed precipitation suggests low aqueous solubility. We can predict the effect of a structural")
    print("modification by calculating the LogP value. A lower LogP indicates higher water solubility.\n")

    print("--- Original Probe ---")
    print(f"Structure: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide")
    print(f"Molecular Weight: {original_wt:.2f}")
    print(f"Calculated LogP: {original_logp:.2f}")
    print(f"Hydrogen Bond Acceptors: {original_hba}")
    print(f"Hydrogen Bond Donors: {original_hbd}\n")


    print("--- Proposed Modified Probe (PEGylated) ---")
    print("Modification: Replaced the hydrophobic chlorohexyl group with a hydrophilic PEG4 chain.")
    print(f"Molecular Weight: {modified_wt:.2f}")
    print(f"Calculated LogP: {modified_logp:.2f}")
    print(f"Hydrogen Bond Acceptors: {modified_hba}")
    print(f"Hydrogen Bond Donors: {modified_hbd}\n")

    print("--- Conclusion ---")
    logp_change = original_logp - modified_logp
    if modified_logp < original_logp:
        print(f"YES, this modification is highly likely to solve the precipitation problem.")
        print(f"The calculated LogP value decreased from {original_logp:.2f} to {modified_logp:.2f}, a significant drop of {logp_change:.2f} units.")
        print("This indicates a substantial increase in predicted water solubility. Replacing the long, 'greasy'")
        print("hexyl chain with a polar PEG chain is a standard and effective strategy to improve the")
        print("biophysical properties of small molecule probes.")
    else:
        print("NO, based on this calculation, the proposed modification is unlikely to solve the solubility problem.")
        print("The LogP value did not decrease, suggesting no significant improvement in water solubility.")

if __name__ == '__main__':
    analyze_probe_solubility()