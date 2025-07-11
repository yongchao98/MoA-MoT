# To run this code, you need to have the RDKit library installed.
# You can install it via conda: conda install -c conda-forge rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it to run this code, e.g., 'conda install -c conda-forge rdkit'")
    exit()

def solve_chemistry_problem():
    """
    This function identifies Compound C and prints its chemical properties.
    """
    # Step 1: Define the structure of Compound C in SMILES format.
    # The structure is deduced as 3-(diethylamino)-1,6,8-trihydroxyxanthylium cation.
    # We use the resonance form with a positive charge on the oxygen atom (oxonium ion).
    smiles_C = "CCN(CC)C1=CC2=C(C=C1O)C(=[O+])C3=C(C2)C=C(O)C=C3O"

    # Step 2: Create a molecule object from the SMILES string.
    mol_C = Chem.MolFromSmiles(smiles_C)

    if mol_C is None:
        print("Could not parse the SMILES string for Compound C.")
        return

    # Add explicit hydrogens to get the correct formula and weight for the cation.
    mol_C_with_Hs = Chem.AddHs(mol_C)
    
    # Step 3: Calculate properties of Compound C.
    compound_name = "3-(diethylamino)-1,6,8-trihydroxyxanthylium"
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol_C_with_Hs)
    # For the cation [C17H18NO4]+
    exact_mass = Descriptors.ExactMolWt(mol_C_with_Hs)

    # Step 4: Print the results.
    print("--- Identification of Compound C ---")
    print(f"Compound Name: {compound_name}")
    print(f"SMILES String: {smiles_C}")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Exact Mass (Monoisotopic): {exact_mass:.4f} Da")
    print("")

    # As requested, outputting all numbers from the reaction description.
    print("--- Numerical Values from Reaction Description ---")
    numbers = [1, 3, 5, 1.04, 70, 2, 0.3, 3, 9, 10, 170, 4]
    print("The numbers mentioned in the reaction steps are:")
    # Print each number as requested
    for num in numbers:
        print(num)

if __name__ == "__main__":
    solve_chemistry_problem()
