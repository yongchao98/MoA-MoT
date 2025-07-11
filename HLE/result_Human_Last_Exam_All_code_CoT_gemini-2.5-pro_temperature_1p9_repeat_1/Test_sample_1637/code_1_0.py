# First, ensure the rdkit library is installed. If not, run:
# pip install rdkit-pypi

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("RDKit library not found. Please install it first by running:")
    print("pip install rdkit-pypi")
    exit()

# Define the chemical structures using their SMILES representation.
# Original Probe: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
original_smiles = "COc1ccc2c(c1)Sc1cc(OCC(=O)NCCOCCOCCCCCCCl)ccc1C2=O"

# Modified Probe: Here we replace the -C(=O)NH- amide group with a -CH2-O- ether group
# to make the linker more flexible and PEG-like.
modified_smiles = "COc1ccc2c(c1)Sc1cc(OCCOCCOCCOCCCCCCCl)ccc1C2=O"

# Create molecule objects from the SMILES strings
original_mol = Chem.MolFromSmiles(original_smiles)
modified_mol = Chem.MolFromSmiles(modified_smiles)

if original_mol is None or modified_mol is None:
    print("Error: Could not generate a molecule from one of the provided SMILES strings.")
else:
    # --- Calculations for the Original Probe ---
    original_mw = Descriptors.MolWt(original_mol)
    original_logp = Descriptors.MolLogP(original_mol)
    original_formula = rdMolDescriptors.CalcMolFormula(original_mol)

    # --- Calculations for the Modified Probe ---
    modified_mw = Descriptors.MolWt(modified_mol)
    modified_logp = Descriptors.MolLogP(modified_mol)
    modified_formula = rdMolDescriptors.CalcMolFormula(modified_mol)

    # --- Print the comparison and final results ---
    print("Analysis of Proposed Molecular Modification to Improve Solubility")
    print("-" * 65)
    print("Property                | Original Probe          | Modified Probe")
    print("-" * 65)
    print(f"Molecular Formula       | {original_formula:<23} | {modified_formula}")
    print(f"Molecular Weight        | {original_mw:<23.2f} | {modified_mw:.2f}")
    print(f"Hydrophobicity (cLogP)  | {original_logp:<23.2f} | {modified_logp:.2f}")
    print("-" * 65)
    
    # Fulfilling the request to output numbers in a final equation
    print("\nEquation for Change in Hydrophobicity (cLogP):")
    print(f"{modified_logp:.2f} (Modified cLogP) - {original_logp:.2f} (Original cLogP) = {modified_logp - original_logp:.2f}")
    
    if modified_logp < original_logp:
        print("\nConclusion: The modification lowers the cLogP value. This indicates a decrease in")
        print("hydrophobicity and a likely increase in water solubility. This change should help")
        print("solve the precipitation problem.")
    else:
        print("\nConclusion: The modification does not significantly lower the cLogP value and is")
        print("unlikely to resolve the precipitation issue.")
