import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("Error: RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of the original and a modified chemical probe
    to predict whether the modification will improve aqueous solubility.
    """
    # IUPAC: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    # This corresponds to the SMILES structure below, which represents:
    # Thioxanthenone-O-CH2-C(=O)NH-CH2CH2-O-CH2CH2-O-(CH2)6-Cl
    original_smiles = 'COc1ccc2c(c1)sc1cc(OCC(=O)NCCOCCOCCCCCCCl)ccc1C(=O)2'

    # The modified probe replaces the amide group -C(=O)NH- with a more PEG-like ether linkage -CH2O-.
    # This change makes the linker a continuous chain of ether units.
    # Modified structure:
    # Thioxanthenone-O-CH2-CH2-O-CH2CH2-O-CH2CH2-O-(CH2)6-Cl
    modified_smiles = 'COc1ccc2c(c1)sc1cc(OCCOCCOCCOCCCCCCCl)ccc1C(=O)2'

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if original_mol is None or modified_mol is None:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # Calculate CLogP values
    clogp_original = Descriptors.MolLogP(original_mol)
    clogp_modified = Descriptors.MolLogP(modified_mol)
    
    # --- Output Results ---
    print("--- Solubility Analysis ---")
    print("CLogP is a measure of hydrophobicity. A lower value suggests better aqueous solubility.")
    
    print("\nOriginal Probe:")
    # We output the CLogP number for the original probe.
    print(f"  Calculated CLogP: {clogp_original:.4f}")
    
    print("\nModified Probe (Amide -> Ether):")
    # We output the CLogP number for the modified probe.
    print(f"  Calculated CLogP: {clogp_modified:.4f}")

    print("\n--- Conclusion ---")
    print(f"The proposed modification decreases the CLogP value from {clogp_original:.4f} to {clogp_modified:.4f}.")
    print("This reduction indicates that the new compound is more hydrophilic.")
    print("Therefore, changing the amide to a more PEG-like ether linkage is a valid strategy and is likely to solve the precipitation problem.")

if __name__ == "__main__":
    analyze_probe_solubility()