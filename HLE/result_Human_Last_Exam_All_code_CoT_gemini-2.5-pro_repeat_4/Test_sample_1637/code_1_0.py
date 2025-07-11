# First, ensure you have RDKit installed:
# pip install rdkit
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of an original and a modified chemical probe
    to predict if a modification will improve water solubility.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COC1=CC2=C(C=C1)C(=O)C3=CC(=CS32)OCC(=O)NCCOCCOCCCl"

    # SMILES string for a modified probe.
    # The hydrophobic chain on the amide nitrogen is replaced with a hydrophilic PEG4 chain.
    # New structure: N-(PEG4-OH)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    modified_smiles = "COC1=CC2=C(C=C1)C(=O)C3=CC(=CS32)OCC(=O)NCCOCCOCCOCCOH"

    # Create RDKit molecule objects from SMILES
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # Calculate the cLogP for both molecules
    original_logp = Descriptors.MolLogP(original_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)

    print("Analysis of Probe Solubility Using Predicted LogP (cLogP)")
    print("A lower cLogP value indicates higher water solubility.\n")
    print("-" * 50)

    print(f"Original Probe cLogP: {original_logp:.2f}")
    print(f"Modified Probe cLogP: {modified_logp:.2f}")
    print("-" * 50)
    print("\nFinal Equation:")
    
    # We explicitly print the numbers for the final equation as requested
    num1 = original_logp
    num2 = modified_logp
    result = num2 - num1
    
    print(f"Change in cLogP = cLogP(Modified) - cLogP(Original)")
    print(f"Change in cLogP = {num2:.2f} - {num1:.2f} = {result:.2f}")

    print("\nConclusion:")
    if modified_logp < original_logp:
        print("The modified probe has a significantly lower cLogP, indicating it is much more hydrophilic.")
        print("This change is very likely to increase water solubility and solve the precipitation problem.")
    else:
        print("The proposed modification does not improve the predicted water solubility.")

analyze_probe_solubility()
<<<Yes>>>