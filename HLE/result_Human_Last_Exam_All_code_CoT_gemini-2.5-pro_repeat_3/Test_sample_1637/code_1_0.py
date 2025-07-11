# You may need to install the RDKit library first:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen

def analyze_probe_solubility():
    """
    Analyzes and compares the hydrophobicity of the original probe
    and a modified, more soluble version.
    """
    # Structure of the original probe:
    # N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
    original_smiles = "COc1ccc2c(c1)Sc3cc(OCC(=O)NCCOCCOCCCCCCCl)ccc3C2=O"

    # A proposed new structure where the hydrophobic -NH-(linker)-Cl part is replaced
    # with a more hydrophilic PEG linker, -NH-(CH2CH2O)4-H. This drastically
    # increases water solubility by replacing the greasy hexyl chain with a PEG4 chain.
    modified_smiles = "COc1ccc2c(c1)Sc3cc(OCC(=O)NCCOCCOCCOCCOH)ccc3C2=O"

    try:
        # Create RDKit molecule objects from SMILES strings
        original_mol = Chem.MolFromSmiles(original_smiles)
        modified_mol = Chem.MolFromSmiles(modified_smiles)

        if original_mol is None or modified_mol is None:
            print("Error: Invalid SMILES string provided. Please check the structures.")
            return

        # Calculate properties: cLogP (hydrophobicity) and Molecular Weight
        original_clogp = Crippen.MolLogP(original_mol)
        modified_clogp = Crippen.MolLogP(modified_mol)

        original_mw = Descriptors.MolWt(original_mol)
        modified_mw = Descriptors.MolWt(modified_mol)

        print("--- Analysis of Probe Solubility ---")
        print("A lower cLogP value suggests higher hydrophilicity and better water solubility.\n")

        print("Original Probe:")
        print(f"  Molecular Weight: {original_mw:.2f}")
        print(f"  cLogP (Hydrophobicity): {original_clogp:.2f}\n")

        print("Proposed PEGylated Probe:")
        print(f"  Molecular Weight: {modified_mw:.2f}")
        print(f"  cLogP (Hydrophobicity): {modified_clogp:.2f}\n")

        print("--- Conclusion ---")
        print("The proposed modification with a PEG linker significantly reduces the cLogP value,")
        print("predicting a substantial increase in aqueous solubility.\n")

        # Fulfilling the request to show the numbers in a final equation
        print("Equation for Change in Hydrophobicity (cLogP):")
        reduction = original_clogp - modified_clogp
        print(f"Original cLogP - Modified cLogP = Total Reduction")
        print(f"{original_clogp:.2f} - ({modified_clogp:.2f}) = {reduction:.2f}")
        print(f"\nThe hydrophobicity is predicted to decrease by a factor of ~10^{int(round(reduction))}, which is a dramatic improvement.")

    except ImportError:
        print("RDKit library not found. Please install it using 'pip install rdkit'")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    analyze_probe_solubility()