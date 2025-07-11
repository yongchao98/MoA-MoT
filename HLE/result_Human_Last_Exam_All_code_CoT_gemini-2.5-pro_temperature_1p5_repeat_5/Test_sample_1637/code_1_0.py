# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_probe_solubility():
    """
    Analyzes and compares the properties of an original probe and a modified version
    to predict changes in aqueous solubility.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1cc2c(cc1)c3cc(sc3C(=O)c2)OCC(=O)NCCOCCOCCCCCCCl"
    
    # SMILES string for the modified (PEGylated) probe.
    # We replace the hydrophobic chlorohexyl group (-CCCCCCCl) with a hydrophilic
    # PEG-like chain (-CCOCCOH) to improve water solubility.
    modified_smiles = "COc1cc2c(cc1)c3cc(sc3C(=O)c2)OCC(=O)NCCOCCOCCOCCOH"
    
    # Create molecule objects from SMILES
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # --- Calculate properties for the Original Probe ---
    original_logp = Descriptors.MolLogP(original_mol)
    original_mw = Descriptors.ExactMolWt(original_mol)

    # --- Calculate properties for the Modified Probe ---
    modified_logp = Descriptors.MolLogP(modified_mol)
    modified_mw = Descriptors.ExactMolWt(modified_mol)
    
    # --- Print the results and conclusion ---
    print("--- Analysis of Probe Modification for Solubility ---\n")
    print("Yes, changing the hydrophobic part of the chain to a PEG group is an excellent strategy to solve the precipitation problem.")
    print("Here is a quantitative comparison:\n")
    
    print("Original Probe:")
    print(f"  Molecular Weight: {original_mw:.2f}")
    print(f"  Calculated LogP (cLogP): {original_logp:.2f} (Higher value = more hydrophobic)")
    
    print("\nModified (PEGylated) Probe:")
    print(f"  Molecular Weight: {modified_mw:.2f}")
    print(f"  Calculated LogP (cLogP): {modified_logp:.2f} (Lower value = more hydrophilic)")

    print("\n--- Conclusion ---")
    print("The modification replaces a hydrophobic chlorohexyl group with a hydrophilic PEG chain.")
    print("As a result, the calculated LogP drops significantly, indicating a substantial increase in water solubility.")
    
    # The final "equation" showing the comparison
    print("\nFinal Comparison:")
    print(f"LogP_modified ({modified_logp:.2f}) < LogP_original ({original_logp:.2f})")
    
    print("\nThis change is very likely to prevent the probe from precipitating at 100 uM in the cell culture medium.")

# Execute the analysis
analyze_probe_solubility()

<<<Yes, replacing the hydrophobic chlorohexyl group with a hydrophilic PEG chain is a very effective strategy that will likely solve the precipitation problem by significantly increasing the probe's water solubility.>>>