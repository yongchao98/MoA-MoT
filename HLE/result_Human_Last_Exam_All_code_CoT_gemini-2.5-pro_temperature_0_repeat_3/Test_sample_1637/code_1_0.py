# You may need to install the RDKit library first:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Crippen

def analyze_probe_solubility():
    """
    Calculates and compares the LogP values for the original probe and
    a modified version with a PEG chain to estimate the change in solubility.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1cc2c(sc3cc(OCC(=O)NCCOCCOCCCCCCCl)ccc3C2=O)cc1"

    # SMILES string for a proposed modified probe where the entire tail
    # is replaced with a PEG8 (8 ethylene glycol units) chain.
    # This is a common strategy to improve solubility.
    modified_smiles = "COc1cc2c(sc3cc(OCCOCCOCCOCCOCCOCCOCCOCCO)ccc3C2=O)cc1"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings. Please check them.")
        return

    # Calculate the Crippen LogP value for both molecules
    original_logp = Crippen.MolLogP(original_mol)
    modified_logp = Crippen.MolLogP(modified_mol)

    print("Analysis of Probe Solubility (estimated by LogP):")
    print("A lower LogP value suggests higher water solubility.\n")

    # Output the numbers for the "equation": Original vs. Modified
    print(f"Original Probe LogP = {original_logp:.2f}")
    print(f"Modified Probe LogP = {modified_logp:.2f}")
    print("-" * 30)
    print(f"Change in LogP = {modified_logp - original_logp:.2f}\n")

    if modified_logp < original_logp:
        print("Conclusion: The modification significantly decreases the LogP value.")
        print("This strongly suggests the PEGylated probe will have better water solubility,")
        print("which should help solve the precipitation problem.")
    else:
        print("Conclusion: The proposed modification does not appear to improve hydrophilicity.")

# Run the analysis
analyze_probe_solubility()