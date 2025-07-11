import sys
# First, we ensure the necessary rdkit library is available.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it before running this script using the command:")
    print("pip install rdkit-pypi")
    sys.exit(1)

# We define the chemical structures using the SMILES notation.
# Original Probe: N-(2-(2-((6-chlorohexyl)oxy)ethoxy)ethyl)-2-((7-methoxy-9-oxo-9H-thioxanthen-2-yl)oxy)acetamide
original_smiles = 'COC1=CC2=C(C=C1)SC3=CC(OCC(=O)NCCOCCOCCCCCCCl)=CC=C3C2=O'

# Modified Probe with one additional ethylene glycol unit to make it more "PEG-like".
modified_smiles = 'COC1=CC2=C(C=C1)SC3=CC(OCC(=O)NCCOCCOCCOCCCCCCCl)=CC=C3C2=O'

# We create molecule objects from the SMILES strings for calculation.
original_mol = Chem.MolFromSmiles(original_smiles)
modified_mol = Chem.MolFromSmiles(modified_smiles)

if not original_mol or not modified_mol:
    print("Error: Invalid chemical structure format (SMILES).")
    sys.exit(1)

# --- Perform calculations for the original probe ---
original_logp = Descriptors.MolLogP(original_mol)
original_formula = rdMolDescriptors.CalcMolFormula(original_mol)
original_tpsa = Descriptors.TPSA(original_mol)

# --- Perform calculations for the modified probe ---
modified_logp = Descriptors.MolLogP(modified_mol)
modified_formula = rdMolDescriptors.CalcMolFormula(modified_mol)
modified_tpsa = Descriptors.TPSA(modified_mol)

# --- Print the analysis and comparison ---
print("--- Solubility Analysis for the Probe ---")
print("A lower LogP value (hydrophobicity) and higher TPSA (polarity) suggest better water solubility.\n")

print(f"--- Original Probe ---")
print(f"Formula: {original_formula}")
print(f"Calculated LogP: {original_logp:.2f}")
print(f"Topological Polar Surface Area (TPSA): {original_tpsa:.2f} Å²\n")

print(f"--- Modified Probe (with one extra PEG unit) ---")
print(f"Formula: {modified_formula}")
print(f"Calculated LogP: {modified_logp:.2f}")
print(f"Topological Polar Surface Area (TPSA): {modified_tpsa:.2f} Å²\n")

# --- Conclude with a final equation and recommendation ---
print("--- Comparison and Final Equation ---")
logp_change = modified_logp - original_logp
tpsa_change = modified_tpsa - original_tpsa

print("The change in hydrophobicity (LogP) can be expressed with the equation:")
print(f"{modified_logp:.2f} (modified) - {original_logp:.2f} (original) = {logp_change:+.2f}\n")

print("The change in polarity (TPSA) is:")
print(f"{modified_tpsa:.2f} (modified) - {original_tpsa:.2f} (original) = {tpsa_change:+.2f} Å²\n")


if logp_change < 0:
    print("Conclusion: Yes, this modification is likely to solve the solubility problem.")
    print("Adding a PEG unit decreases the LogP and increases the TPSA, both of which strongly suggest improved water solubility.")
else:
    print("Conclusion: No, this modification is unlikely to be sufficient to solve the solubility problem.")

<<<Yes>>>