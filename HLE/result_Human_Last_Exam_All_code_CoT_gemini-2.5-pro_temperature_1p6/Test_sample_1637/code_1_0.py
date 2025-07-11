# The rdkit library is required to run this script.
# If you don't have it, you can install it by running: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("Error: The 'rdkit' library is not installed.")
    print("Please install it using the command: pip install rdkit")
    exit()

def analyze_probe_solubility():
    """
    Calculates and compares the LogP values of an original and modified chemical probe
    to predict the effect of a structural change on aqueous solubility.
    """

    # Step 1: Define the chemical structures using SMILES strings.
    # Original Probe: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    # Structure: [Thioxanthenone]-O-CH2-C(=O)-NH-[PEG-Cl_chain]
    original_smiles = "COc1ccc2c(c1)sc1ccc(OCC(=O)NCCOCCOCCCCCCCl)cc1C2=O"
    original_mol = Chem.MolFromSmiles(original_smiles)

    # Proposed Modified Probe: Replaces the amide linker (-C(=O)NH-) with a more flexible and
    # hydrophilic ether linker, effectively extending the PEG chain.
    # Structure: [Thioxanthenone]-O-CH2-CH2-O-[Longer_PEG-Cl_chain]
    modified_smiles = "COc1ccc2c(c1)sc1ccc(OCCOCCOCCOCCCCCCCl)cc1C2=O"
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if original_mol is None or modified_mol is None:
        print("Error: Failed to create molecule from SMILES string. Check definitions.")
        return

    # Step 2: Calculate the partition coefficient (LogP) for both molecules.
    # Crippen's LogP is a reliable method for this estimation.
    original_logp = Descriptors.MolLogP(original_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)

    # Step 3: Print the results and provide a conclusion.
    print("--- Analysis of Probe Solubility ---")
    print("A lower LogP value suggests higher hydrophilicity and better solubility in aqueous media.\n")

    print("Original Probe:")
    print(f"  SMILES: {original_smiles}")
    print(f"  Calculated LogP: {original_logp:.2f}\n")

    print("Modified Probe (with PEG-like ether linker):")
    print(f"  SMILES: {modified_smiles}")
    print(f"  Calculated LogP: {modified_logp:.2f}\n")

    print("--- Conclusion ---")
    if modified_logp < original_logp:
        print(f"The modified probe's LogP of {modified_logp:.2f} is lower than the original's LogP of {original_logp:.2f}.")
        print("This calculation strongly suggests that replacing the amide with a PEG-like ether linker will increase the probe's water solubility.")
        print("\nTherefore, this modification is a promising strategy to solve the precipitation problem.")
    else:
        print(f"The modified probe's LogP ({modified_logp:.2f}) is not lower than the original's ({original_logp:.2f}).")
        print("Based on this calculation, this specific change may not solve the solubility issue.")

    print("\nDisclaimer: This is a computational prediction. Experimental validation is required to confirm the change in solubility and to ensure the probe's function is not negatively affected.")

if __name__ == '__main__':
    analyze_probe_solubility()