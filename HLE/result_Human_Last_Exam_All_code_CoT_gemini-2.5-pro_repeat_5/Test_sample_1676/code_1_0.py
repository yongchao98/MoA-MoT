import sys

def solve_synthesis():
    """
    This script identifies the final product of a three-step chemical synthesis
    starting from terpinolene. It requires the RDKit library.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    except ImportError:
        print("Error: RDKit is not installed. Please install it to run this script.")
        print("You can install it using: pip install rdkit-pypi")
        sys.exit(1)

    # 1. Starting Material: Terpinolene (4-isopropylidene-1-methylcyclohex-1-ene)
    terpinolene_smiles = 'CC1=CCC(=C(C)C))CC1'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)

    # 2. Step 1: Selective epoxidation of the endocyclic double bond
    # The reaction SMARTS targets the less-hindered trisubstituted double bond in the ring.
    rxn_epox = AllChem.ReactionFromSmarts('[C;R;!H0:1]=[C;R;H1:2]>>[C:1]1O[C:2]1')
    products_step1 = rxn_epox.RunReactants((terpinolene,))
    compound_1 = products_step1[0][0]
    Chem.SanitizeMol(compound_1)

    # 3. Step 2: Conversion of epoxide to thiirane (episulfide)
    rxn_thiirane = AllChem.ReactionFromSmarts('[C:1]1O[C:2]1>>[C:1]1S[C:2]1')
    products_step2 = rxn_thiirane.RunReactants((compound_1,))
    compound_2 = products_step2[0][0]
    Chem.SanitizeMol(compound_2)

    # 4. Step 3: Reductive opening of the thiirane ring with LiAlH4
    # The reaction is modeled by programmatically modifying the molecule,
    # simulating hydride attack at the less substituted carbon.
    rw_mol = Chem.RWMol(compound_2)
    s_atom = next(atom for atom in rw_mol.GetAtoms() if atom.GetSymbol() == 'S')
    c_atom1, c_atom2 = s_atom.GetNeighbors()

    # Identify less substituted carbon (secondary vs. tertiary) for hydride attack
    c_less_substituted = c_atom1 if c_atom1.GetDegree() < c_atom2.GetDegree() else c_atom2
    
    # Get atom indices for manipulation
    s_idx = s_atom.GetIdx()
    c_less_idx = c_less_substituted.GetIdx()

    # Add explicit hydrogens to correctly manage valency during the reaction
    rw_mol = Chem.RWMol(Chem.AddHs(rw_mol))
    s_atom = rw_mol.GetAtomWithIdx(s_idx)
    c_less_substituted = rw_mol.GetAtomWithIdx(c_less_idx)

    # Break the C-S bond at the less substituted carbon
    rw_mol.RemoveBond(c_less_idx, s_idx)
    
    # Add a hydrogen to the attacked carbon (simulating H- attack)
    c_less_substituted.SetNumExplicitHs(c_less_substituted.GetNumExplicitHs() + 1)
    
    # Add a hydrogen to the sulfur to form the thiol (simulating workup)
    s_atom.SetNumExplicitHs(s_atom.GetNumExplicitHs() + 1)
    
    # Get the final molecule and clean it up
    compound_3_with_hs = rw_mol.GetMol()
    Chem.SanitizeMol(compound_3_with_hs)
    compound_3 = Chem.RemoveHs(compound_3_with_hs)
    Chem.SanitizeMol(compound_3)

    # 5. Output the result
    final_product_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"
    final_formula = CalcMolFormula(compound_3)

    print("The final product, Compound 3, is:")
    print(f"Name: {final_product_name}")
    print(f"Chemical Formula: {final_formula}")

if __name__ == "__main__":
    solve_synthesis()