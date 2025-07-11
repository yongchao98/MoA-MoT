# To run this script, you first need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_compound_info(name, smiles):
    """Generates a dictionary with information about a chemical compound."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"name": name, "error": "Invalid SMILES string"}
    
    formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    exact_mass = Descriptors.ExactMolWt(mol)
    
    return {
        "name": name,
        "smiles": smiles,
        "formula": formula,
        "mass": f"{exact_mass:.2f}"
    }

def solve_chemistry_problem():
    """
    Analyzes the described reaction sequence and identifies Compound 3.
    """
    # Step 1: Define the SMILES strings for each compound in the reaction path.
    # Terpinolene has a trisubstituted endocyclic and a tetrasubstituted exocyclic double bond.
    # m-CPBA epoxidizes the more substituted (exocyclic) double bond.
    smiles_start = "CC1=CCC(=C(C)C)CC1"  # Terpinolene (Starting Material)
    
    # Compound 1 is the epoxide of Terpinolene at the exocyclic double bond.
    smiles_c1 = "CC1=CCC(C2(C)OC2C)CC1"
    
    # Compound 2 is the thiirane (episulfide) analog, from reaction with thioformamide.
    smiles_c2 = "CC1=CCC(C2(C)SC2C)CC1"
    
    # Compound 3 results from the LiAlH4 reduction of the thiirane, which regenerates the alkene.
    # Therefore, Compound 3 is the same as the starting material.
    smiles_c3 = smiles_start
    
    # Step 2: Get information for each compound
    terpinolene = get_compound_info("Terpinolene", smiles_start)
    compound_1 = get_compound_info("Compound 1", smiles_c1)
    compound_2 = get_compound_info("Compound 2", smiles_c2)
    compound_3 = get_compound_info("Compound 3", smiles_c3)
    
    # Step 3: Print the results, including the final "equation" with "numbers" (masses).
    print("The reaction pathway is as follows:")
    print("-" * 30)
    print(f"Starting Material: {terpinolene['name']}")
    print(f"Structure (SMILES): {terpinolene['smiles']}")
    print(f"Formula: {terpinolene['formula']}, Exact Mass: {terpinolene['mass']}")
    print("-" * 30)
    
    print("Step 1: Epoxidation with m-CPBA")
    print(f"Product: {compound_1['name']}")
    print(f"Structure (SMILES): {compound_1['smiles']}")
    print(f"Formula: {compound_1['formula']}, Exact Mass: {compound_1['mass']}")
    print("-" * 30)
    
    print("Step 2: Conversion to Thiirane")
    print(f"Product: {compound_2['name']}")
    print(f"Structure (SMILES): {compound_2['smiles']}")
    print(f"Formula: {compound_2['formula']}, Exact Mass: {compound_2['mass']}")
    print("-" * 30)

    print("Step 3: Reduction with LiAlH4")
    print(f"Final Product: {compound_3['name']}")
    print(f"Structure (SMILES): {compound_3['smiles']}")
    print(f"Formula: {compound_3['formula']}, Exact Mass: {compound_3['mass']}")
    print("-" * 30)
    
    # Final summary equation with molecular formulas and masses as the "numbers"
    print("\nSummary Equation of Transformations:")
    pathway_str = (
        f"{terpinolene['name']} ({terpinolene['formula']}, {terpinolene['mass']}) -> "
        f"{compound_1['name']} ({compound_1['formula']}, {compound_1['mass']}) -> "
        f"{compound_2['name']} ({compound_2['formula']}, {compound_2['mass']}) -> "
        f"{compound_3['name']} ({compound_3['formula']}, {compound_3['mass']})"
    )
    print(pathway_str)
    
    print(f"\nConclusion: Compound 3 is {compound_3['name']}.")

if __name__ == "__main__":
    solve_chemistry_problem()