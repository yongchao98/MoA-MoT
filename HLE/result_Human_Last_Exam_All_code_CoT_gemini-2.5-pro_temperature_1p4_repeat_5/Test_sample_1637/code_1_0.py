import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit")
    sys.exit(1)

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of an original probe and a
    modified version to predict if the modification will improve solubility.
    """

    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COC1=CC2=C(C=C1)SC3=CC=C(OCC(=O)NCCOCCOCCCCCCCl)C=C3C2=O"

    # SMILES string for the modified probe, where the acetamide linker (-O-CH2-C(=O)-NH-)
    # is replaced with a more hydrophilic PEG-like linker (-O-CH2CH2-O-CH2CH2-).
    modified_smiles = "COC1=CC2=C(C=C1)SC3=CC=C(OCCOCCOCCOCCOCCCCCCCl)C=C3C2=O"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # --- Calculate Physicochemical Properties ---

    # 1. LogP (Hydrophobicity)
    original_logp = Descriptors.MolLogP(original_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)

    # 2. Hydrogen Bond Donors
    original_hbd = Descriptors.NumHDonors(original_mol)
    modified_hbd = Descriptors.NumHDonors(modified_mol)

    # 3. Hydrogen Bond Acceptors
    original_hba = Descriptors.NumHAcceptors(original_mol)
    modified_hba = Descriptors.NumHAcceptors(modified_mol)
    
    # 4. Molecular Weight
    original_mw = Descriptors.MolWt(original_mol)
    modified_mw = Descriptors.MolWt(modified_mol)


    # --- Print Results and Conclusion ---
    print("--- Analysis of Probe Solubility ---\n")
    print("Original Probe:")
    print(f"  - Predicted LogP: {original_logp:.2f}")
    print(f"  - Hydrogen Bond Donors: {original_hbd}")
    print(f"  - Hydrogen Bond Acceptors: {original_hba}")
    print(f"  - Molecular Weight: {original_mw:.2f}\n")


    print("Modified Probe (Amide replaced with PEG linker):")
    print(f"  - Predicted LogP: {modified_logp:.2f}")
    print(f"  - Hydrogen Bond Donors: {modified_hbd}")
    print(f"  - Hydrogen Bond Acceptors: {modified_hba}")
    print(f"  - Molecular Weight: {modified_mw:.2f}\n")

    print("--- Conclusion ---")
    print("A lower LogP value indicates higher water solubility.")
    
    # Final Equation Output
    logp_change = original_logp - modified_logp
    print("The final comparison is based on the change in LogP:")
    print(f"LogP_Improvement = LogP(Original) - LogP(Modified)")
    print(f"LogP_Improvement = {original_logp:.2f} - {modified_logp:.2f} = {logp_change:.2f}")

    if modified_logp < original_logp:
        print("\nThe modification significantly decreases the calculated LogP, indicating a substantial increase in hydrophilicity (water solubility).")
        print("This change is very likely to solve the precipitation problem.")
        global conclusion_answer
        conclusion_answer = "<<<Yes>>>"
    else:
        print("\nThe modification does not decrease the calculated LogP. It is unlikely to solve the precipitation problem.")
        conclusion_answer = "<<<No>>>"

# Execute the analysis
analyze_probe_solubility()
print(conclusion_answer)