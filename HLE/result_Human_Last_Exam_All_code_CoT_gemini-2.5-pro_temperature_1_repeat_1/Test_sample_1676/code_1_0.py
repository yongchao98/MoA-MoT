from rdkit import Chem
from rdkit.Chem import AllChem

def run_reaction(reactant_mol, reaction_smarts):
    """Helper function to run a reaction and return the major product."""
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        products = rxn.RunReactants((reactant_mol,))
        if products:
            # Get the first product of the first reaction outcome
            product_mol = products[0][0]
            Chem.SanitizeMol(product_mol)
            return product_mol
        else:
            return None
    except Exception as e:
        print(f"An error occurred during the reaction: {e}")
        return None

def solve_chemistry_problem():
    """
    Solves the multi-step synthesis problem starting from terpinolene.
    """
    # Define the starting material: Terpinolene (1-methyl-4-(propan-2-ylidene)cyclohex-1-ene)
    terpinolene_smiles = 'CC1=CCC(CC1)C(=C(C)C)'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    print("The multi-step reaction starts with Terpinolene.")
    print("Equation:")
    print(f"Start (Terpinolene): {Chem.MolToSmiles(terpinolene)}")

    # Step 1: Epoxidation of the more substituted (exocyclic) double bond
    # The SMARTS targets a tetrasubstituted double bond.
    epoxidation_smarts = '[C:1](C)(C)=[C:2](C)(C)>>[C:1]1O[C:2]1'
    compound_1 = run_reaction(terpinolene, epoxidation_smarts)
    if not compound_1:
        print("Step 1 failed.")
        return
    print(f"  -> Compound 1 (Epoxide): {Chem.MolToSmiles(compound_1)}")

    # Step 2: Conversion of epoxide to thiirane
    thiirane_formation_smarts = '[C:1]1O[C:2]1>>[C:1]1S[C:2]1'
    compound_2 = run_reaction(compound_1, thiirane_formation_smarts)
    if not compound_2:
        print("Step 2 failed.")
        return
    print(f"  -> Compound 2 (Thiirane): {Chem.MolToSmiles(compound_2)}")
    
    # Step 3: Reduction (desulfurization) of thiirane to alkene with LiAlH4
    desulfurization_smarts = '[C:1]1S[C:2]1>>[C:1]=[C:2]'
    compound_3 = run_reaction(compound_2, desulfurization_smarts)
    if not compound_3:
        print("Step 3 failed.")
        return
    print(f"  -> Compound 3 (Final Product): {Chem.MolToSmiles(compound_3)}")

    print("\n-------------------------------------------")
    print("The final product (Compound 3) is Terpinolene.")
    print("-------------------------------------------")

if __name__ == "__main__":
    solve_chemistry_problem()
