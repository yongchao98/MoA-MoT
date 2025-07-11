from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_probe_solubility():
    """
    Analyzes and compares the physicochemical properties of an original and
    a modified chemical probe to predict changes in solubility.
    """

    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1cc2c(sc3cc(OC(=O)CNCCOCCOCCCCCCCl)ccc3C2=O)cc1"
    original_mol = Chem.MolFromSmiles(original_smiles)

    # SMILES string for a modified probe.
    # The chloro-hexyl tail has been replaced with a longer, more hydrophilic PEG chain
    # ending in a hydroxyl group to improve solubility.
    modified_smiles = "COc1cc2c(sc3cc(OC(=O)CNCCOCCOCCOCCOCCOH)ccc3C2=O)cc1"
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one or both of the SMILES strings.")
        return

    # --- Calculate properties for the Original Probe ---
    original_mw = Descriptors.MolWt(original_mol)
    original_logp = Descriptors.MolLogP(original_mol)

    # --- Calculate properties for the Modified Probe ---
    modified_mw = Descriptors.MolWt(modified_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)
    
    # --- Print the results and conclusion ---
    print("Analysis of Probe Properties:")
    print("-" * 30)
    
    # The user requested to output each number in the final equation.
    # Here we present the numbers clearly for comparison.
    print("Original Probe:")
    print(f"Molecular Weight = {original_mw:.2f}")
    print(f"Calculated LogP = {original_logp:.2f}\n")
    
    print("Modified Probe (with enhanced PEG linker):")
    print(f"Molecular Weight = {modified_mw:.2f}")
    print(f"Calculated LogP = {modified_logp:.2f}\n")
    
    print("-" * 30)
    print("Conclusion:")
    logp_change = original_logp - modified_logp
    print(f"The modified probe has a Calculated LogP that is {logp_change:.2f} units lower.")
    print("This indicates it is significantly more hydrophilic (water-loving).")
    print("This change is very likely to solve the precipitation problem.")

# Run the analysis
analyze_probe_solubility()
<<<Yes>>>