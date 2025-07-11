# This script uses the RDKit library to verify the molecular formulas
# of the proposed products against the experimental HRMS data.
# To run this, you may need to install RDKit:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def analyze_product(name, smiles, expected_formula):
    """
    Calculates the molecular formula and mass for a given SMILES string
    and prints a comparison with the expected values.
    """
    print(f"--- Analysis for {name} ---")
    try:
        # Create a molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Could not parse SMILES for {name}")
            return

        # Calculate molecular formula and exact molecular weight
        formula = CalcMolFormula(mol)
        exact_mass = Descriptors.ExactMolWt(mol)

        print(f"Proposed Structure (SMILES): {smiles}")
        print(f"Calculated Molecular Formula: {formula}")
        print(f"Expected Molecular Formula (from HRMS): {expected_formula}")

        if formula == expected_formula:
            print("Result: Calculated formula matches the experimental data.")
        else:
            print("Result: Calculated formula does NOT match the experimental data.")
        print(f"Calculated Exact Mass: {exact_mass:.4f}\n")

    except ImportError:
        print("RDKit library not found. Please install it using 'pip install rdkit'")
    except Exception as e:
        print(f"An error occurred: {e}")

# Define the SMILES strings for the products
# Product A is a 1:1 mixture; both regioisomers have the same formula.
# We use methyl 1,2,5-trimethyl-1H-pyrrole-4-carboxylate as a representative.
smiles_A = "CN1C(=C(C(=C1C)C)C(=O)OC)"
formula_A = "C9H13NO2"

# Product B is the regioselective product from the proline-derived münchnone.
# This is methyl 7-methyl-2,3-dihydro-1H-pyrrolizine-6-carboxylate.
smiles_B = "CC1=C(C(=O)OC)C=C2N1CCC2"
formula_B = "C10H13NO2"

# Run the analysis
analyze_product("Product A", smiles_A, formula_A)
analyze_product("Product B", smiles_B, formula_B)