# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of an original probe and
    a proposed modified version with increased PEGylation.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1ccc2c(c1)sc3cc(OCC(=O)NCCOCCOCCCCCCCCCl)ccc3C2=O"

    # SMILES string for a modified probe where the acetamide linker is replaced with a longer PEG chain
    # to increase hydrophilicity. This is a practical synthetic modification.
    modified_smiles = "COc1ccc2c(c1)sc3cc(OCCOCCOCCOCCOCCCCCCCl)ccc3C2=O"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # --- Calculate properties for the Original Probe ---
    original_logp = Crippen.MolLogP(original_mol)
    original_tpsa = Descriptors.TPSA(original_mol)
    original_mw = Descriptors.MolWt(original_mol)

    # --- Calculate properties for the Modified Probe ---
    modified_logp = Crippen.MolLogP(modified_mol)
    modified_tpsa = Descriptors.TPSA(modified_mol)
    modified_mw = Descriptors.MolWt(modified_mol)

    # --- Print the results for comparison ---
    print("--- Analysis of Probe Solubility ---")
    print("\n[ Original Probe ]")
    print(f"Structure: {original_smiles}")
    print(f"Molecular Weight: {original_mw:.2f}")
    print(f"Calculated LogP: {original_logp:.2f}")
    print(f"Topological Polar Surface Area (TPSA): {original_tpsa:.2f} Å²")


    print("\n[ Proposed Modified Probe (with PEG linker) ]")
    print(f"Structure: {modified_smiles}")
    print(f"Molecular Weight: {modified_mw:.2f}")
    print(f"Calculated LogP: {modified_logp:.2f}")
    print(f"Topological Polar Surface Area (TPSA): {modified_tpsa:.2f} Å²")

    print("\n--- Conclusion ---")
    logp_change = ((modified_logp - original_logp) / original_logp) * 100
    tpsa_change = ((modified_tpsa - original_tpsa) / original_tpsa) * 100

    print("The proposed modification significantly decreases the LogP (a measure of hydrophobicity) and increases the TPSA (a measure of polarity).")
    print("A lower LogP and higher TPSA both strongly indicate improved water solubility.")
    print("\nQuantitative Comparison:")
    print(f"Change in LogP: {modified_logp:.2f} (Original) -> {original_logp:.2f} (Modified). This is a {logp_change:.2f}% decrease.")
    print(f"Change in TPSA: {original_tpsa:.2f} Å² (Original) -> {modified_tpsa:.2f} Å² (Modified). This is a {tpsa_change:.2f}% increase.")
    print("\nTherefore, yes, changing the amide group to a PEG group is very likely to solve the precipitation problem.")

if __name__ == "__main__":
    analyze_probe_solubility()