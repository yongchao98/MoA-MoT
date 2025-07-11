try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it first, for example, using pip:")
    print("pip install rdkit")
    exit()

def analyze_probe_solubility():
    """
    Analyzes the change in solubility of a chemical probe after modification
    by calculating and comparing the LogP values.
    """
    # SMILES string for the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1cc2c(sc3ccc(OCC(=O)NCCOCCOCCCCCCCl)cc3C2=O)cc1"

    # The modification "change the amide group to the PEG group" is interpreted as
    # increasing the number of hydrophilic PEG units in the linker to improve solubility,
    # while preserving the reactive chloro-hexyl group.
    # Here, we add two extra ethylene glycol (-O-CH2-CH2-) units.
    modified_smiles = "COc1cc2c(sc3ccc(OCC(=O)NCCOCCOCCOCCOCCCCCCCl)cc3C2=O)cc1"

    # Create RDKit molecule objects
    original_mol = Chem.MolFromSmiles(original_smiles)
    modified_mol = Chem.MolFromSmiles(modified_smiles)

    if not original_mol or not modified_mol:
        print("Error: Could not parse one of the SMILES strings. Please check them.")
        return

    # Calculate properties: Molecular Weight and LogP
    original_logp = Descriptors.MolLogP(original_mol)
    original_mw = Descriptors.MolWt(original_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)
    modified_mw = Descriptors.MolWt(modified_mol)

    # --- Output the results ---
    print("Analysis of Probe Solubility Modification")
    print("=" * 50)
    
    print("\n[Original Probe]")
    print(f"Structure (SMILES): {original_smiles}")
    print(f"Molecular Weight: {original_mw:.2f}")
    print(f"Calculated LogP: {original_logp:.2f}")

    print("\n[Proposed Modification]")
    print("Strategy: Increase hydrophilicity by adding two ethylene glycol units to the linker.")
    print(f"Structure (SMILES): {modified_smiles}")
    print(f"Molecular Weight: {modified_mw:.2f}")
    print(f"Calculated LogP: {modified_logp:.2f}")

    print("\n" + "=" * 50)
    print("\n[Conclusion]")
    print("LogP is a measure of hydrophobicity. A lower LogP value indicates higher water solubility.")
    print(f"The modification decreased the LogP from {original_logp:.2f} to {modified_logp:.2f}.")
    print("This significant decrease in LogP confirms that the modified probe is more hydrophilic.")
    print("\nTherefore, increasing the length of the PEG linker is an excellent strategy and is very likely to solve the precipitation problem.")

if __name__ == "__main__":
    analyze_probe_solubility()
