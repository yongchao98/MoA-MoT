# The user needs to have rdkit installed. It can be installed via conda:
# conda install -c conda-forge rdkit

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem import AllChem

def get_molecule_info(name, smiles_string):
    """Helper function to create a molecule and get its info."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return f"Could not parse SMILES for {name}", None, None
    
    # Add explicit hydrogens for accurate mass calculation and viewing
    mol = Chem.AddHs(mol)
    
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Note: RDKit does not have a built-in IUPAC name generator.
    # The name is provided from chemical knowledge.
    return mol, formula

def solve_chemistry_problem():
    """
    Solves the multi-step synthesis problem and prints the results.
    """
    # Step 0: Starting Material - Terpinolene
    terpinolene_smiles = "CC1=CCC(=C(C)C)CC1"
    terpinolene_mol, terpinolene_formula = get_molecule_info("Terpinolene", terpinolene_smiles)

    # Step 1: Epoxidation of Terpinolene -> Compound 1
    # The more substituted exocyclic double bond is epoxidized.
    compound1_smiles = "CC1=CCC(C2(C)C)O2)CC1"
    compound1_mol, compound1_formula = get_molecule_info("Compound 1", compound1_smiles)

    # Step 2: Epoxide to Thiirane -> Compound 2
    # The epoxide is converted to a thiirane.
    compound2_smiles = "CC1=CCC(C2(C)C)S2)CC1"
    compound2_mol, compound2_formula = get_molecule_info("Compound 2", compound2_smiles)

    # Step 3: Reduction of Thiirane -> Compound 3
    # The thiirane is reductively opened to a thiol.
    compound3_smiles = "CC1=CCC(C(C)(C)S)CC1"
    compound3_mol, compound3_formula = get_molecule_info("Compound 3", compound3_smiles)

    # --- Output the results ---
    print("--- Reaction Sequence ---")
    print(f"Step 0: Starting Material: Terpinolene")
    print(f"   - Molecular Formula: {terpinolene_formula}")
    print(f"   - SMILES: {terpinolene_smiles}\n")

    print(f"Step 1: Epoxidation -> Compound 1")
    print(f"   - Molecular Formula: {compound1_formula}")
    print(f"   - SMILES: {compound1_smiles}\n")

    print(f"Step 2: Thiirane Formation -> Compound 2")
    print(f"   - Molecular Formula: {compound2_formula}")
    print(f"   - SMILES: {compound2_smiles}\n")

    print(f"Step 3: Reduction -> Compound 3 (Final Product)")
    print(f"   - Molecular Formula: {compound3_formula}")
    print(f"   - SMILES: {compound3_smiles}\n")
    
    print("--- Final Product Details (Compound 3) ---")
    iupac_name = "4-(2-mercaptopropan-2-yl)-1-methylcyclohex-1-ene"
    common_name = "p-Menth-1-en-8-thiol"
    print(f"IUPAC Name: {iupac_name}")
    print(f"Common Name: {common_name}")
    
    print("\n--- Final Equation Numbers ---")
    # As requested, printing each number from the final molecular formula.
    # Final formula is C10H18S
    print("The final molecular formula is C10H18S.")
    print("Carbon atoms: 10")
    print("Hydrogen atoms: 18")
    print("Sulfur atoms: 1")

if __name__ == '__main__':
    solve_chemistry_problem()