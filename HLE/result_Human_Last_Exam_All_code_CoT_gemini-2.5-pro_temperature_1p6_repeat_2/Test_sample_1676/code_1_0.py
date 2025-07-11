import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def describe_compound(name, smiles):
    """Creates an RDKit molecule and returns its description."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"{name}\n  SMILES: {smiles}\n  Error: Invalid SMILES string"
    
    formula = CalcMolFormula(mol)
    mol_weight = Descriptors.MolWt(mol)
    
    return (f"{name}\n"
            f"  SMILES: {smiles}\n"
            f"  Molecular Formula: {formula}\n"
            f"  Molecular Weight: {mol_weight:.2f}")

# Define the compounds in the reaction sequence
# The user wants each numbered compound in the final equation.
compounds = {
    "Start (Terpinolene)": "CC1=CCC(=C(C)C)CC1",
    "Compound 1 (Epoxide)": "C1(C)C2CCC(=C(C)C)CC2O1",
    "Compound 2 (Thiirane)": "C1(C)C2CCC(=C(C)C)CC2S1",
    "Compound 3 (Final Product)": "CC1(S)CCC(=C(C)C)CC1"
}

# --- Main Execution ---
print("--- Reaction Pathway Analysis ---")
print("\nStep 1: Epoxidation of Terpinolene")
print(describe_compound("Terpinolene", compounds["Start (Terpinolene)"]))
print("     |")
print("     + m-CPBA, DCM, 0 C")
print("     V")
print(describe_compound("Compound 1", compounds["Compound 1 (Epoxide)"]))

print("\nStep 2: Conversion to Thiirane")
print("Compound 1")
print("     |")
print("     + N,N-dimethylthioformamide, TFA, 60 C")
print("     V")
print(describe_compound("Compound 2", compounds["Compound 2 (Thiirane)"]))

print("\nStep 3: Reduction of Thiirane")
print("Compound 2")
print("     |")
print("     + LiAlH4, 0 C")
print("     V")
print("Compound 3")

print("\n--- Final Product Details ---")
final_product_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"
common_name = "p-menth-8-en-1-thiol"
final_product_smiles = compounds["Compound 3 (Final Product)"]

mol_final = Chem.MolFromSmiles(final_product_smiles)
final_formula = CalcMolFormula(mol_final)
final_mol_weight = Descriptors.MolWt(mol_final)

print(f"The final product, Compound 3, has been identified.\n")
print(f"Systematic Name: {final_product_name}")
print(f"Common Name: {common_name}")
print(f"SMILES: {final_product_smiles}")
print(f"Molecular Formula: {final_formula}")
print(f"Molecular Weight: {final_mol_weight:.2f}")
