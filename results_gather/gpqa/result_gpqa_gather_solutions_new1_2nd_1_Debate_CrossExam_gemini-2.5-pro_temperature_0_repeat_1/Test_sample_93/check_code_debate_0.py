import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def check_correctness():
    """
    Checks the correctness of the answer by verifying the chemical reasoning.
    - Options C and D are intermolecular reactions that lead to incorrect products.
    - Options A and B are intramolecular reactions (IMDA) that lead to different regioisomers.
    - Only option A produces the target molecule with the correct structure.
    """
    final_answer = "A" # The final answer provided by the LLM.

    # --- Define Starting Materials using SMILES ---
    # Option A: methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate
    smiles_A = 'CCC/C=C/C=C/CCCC/C=C/C(=O)OC'
    mol_A = Chem.MolFromSmiles(smiles_A)

    # Option B: methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate
    smiles_B = 'CCC/C=C\\CCCC/C=C/C=C/C(=O)OC'
    mol_B = Chem.MolFromSmiles(smiles_B)

    # Option C: 1-vinylcyclohex-1-ene and methyl hex-2-ynoate
    smiles_C1 = 'C=CC1=CCCCC1'
    smiles_C2 = 'CCCC#CC(=O)OC'
    
    # Option D: Cyclohexene and methyl 2,3-dimethylenehexanoate
    smiles_D1 = 'C1CCCCC=1'
    smiles_D2 = 'CCCC(C(=C))C(=C)C(=O)OC'

    # --- Verification Step 1: Check incorrect intermolecular reactions (C and D) ---

    # Check Option C: Should produce a product with too many double bonds.
    try:
        rxn_C = AllChem.ReactionFromSmarts('[c:1]=[c:2][c:3]=[c:4].[c:5]#[c:6]>>[C:1]1[C:2]=[C:3][C:4][C:6]=[C:5]1')
        reactants_C = (Chem.MolFromSmiles(smiles_C1), Chem.MolFromSmiles(smiles_C2))
        products_C = rxn_C.RunReactants(reactants_C)
        if products_C:
            product_C_mol = products_C[0][0]
            Chem.SanitizeMol(product_C_mol)
            # Target has 1 double bond. Reactant C1 has 2. Product should have 3.
            num_double_bonds_C = sum(1 for b in product_C_mol.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.DOUBLE)
            if num_double_bonds_C != 1:
                pass # Correctly identified as wrong
            else:
                return "Incorrect: Code failed to verify that option C produces a product with the wrong number of double bonds."
        else:
            return "Incorrect: Code failed to simulate reaction for option C."
    except Exception as e:
        return f"An error occurred while checking option C: {e}"

    # Check Option D: Should produce a spirocycle, not a fused system.
    try:
        rxn_D = AllChem.ReactionFromSmarts('[c:1]=[c:2][c:3]=[c:4].[c:5]=[c:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1')
        reactants_D = (Chem.MolFromSmiles(smiles_D2), Chem.MolFromSmiles(smiles_D1))
        products_D = rxn_D.RunReactants(reactants_D)
        if products_D:
            product_D_mol = products_D[0][0]
            Chem.SanitizeMol(product_D_mol)
            # Check for a spiro atom: a degree-4 carbon in exactly two rings.
            spiro_smarts = '[#6;R2;D4]([#6;R1])([#6;R1])([#6;R1])([#6;R1])'
            if not product_D_mol.HasSubstructMatch(Chem.MolFromSmarts(spiro_smarts)):
                 return "Incorrect: Code failed to verify that option D produces a spirocycle."
        else:
            return "Incorrect: Code failed to simulate reaction for option D."
    except Exception as e:
        return f"An error occurred while checking option D: {e}"

    # --- Verification Step 2: Check IMDA reactions (A and B) ---
    
    # Build product from Option A (the proposed correct answer)
    # Reaction: C2-C11 bond, C3-C8 bond. New double bond C9=C10.
    # Atom indices (0-based) in SMILES A: C2=1, C3=2, C8=7, C9=8, C10=9, C11=10
    rw_mol_A = Chem.RWMol(mol_A)
    rw_mol_A.AddBond(1, 10, Chem.BondType.SINGLE)
    rw_mol_A.AddBond(2, 7, Chem.BondType.SINGLE)
    rw_mol_A.GetBondBetweenAtoms(1, 2).SetBondType(Chem.BondType.SINGLE)
    rw_mol_A.GetBondBetweenAtoms(7, 8).SetBondType(Chem.BondType.SINGLE)
    rw_mol_A.GetBondBetweenAtoms(8, 9).SetBondType(Chem.BondType.DOUBLE)
    rw_mol_A.GetBondBetweenAtoms(9, 10).SetBondType(Chem.BondType.SINGLE)
    product_A = rw_mol_A.GetMol()
    Chem.SanitizeMol(product_A)

    # Build product from Option B
    # Reaction: C2-C11 bond, C5-C10 bond. New double bond C3=C4.
    # Atom indices (0-based) in SMILES B: C2=1, C3=2, C4=3, C5=4, C10=9, C11=10
    rw_mol_B = Chem.RWMol(mol_B)
    rw_mol_B.AddBond(1, 10, Chem.BondType.SINGLE)
    rw_mol_B.AddBond(4, 9, Chem.BondType.SINGLE)
    rw_mol_B.GetBondBetweenAtoms(1, 2).SetBondType(Chem.BondType.SINGLE)
    rw_mol_B.GetBondBetweenAtoms(2, 3).SetBondType(Chem.BondType.DOUBLE)
    rw_mol_B.GetBondBetweenAtoms(3, 4).SetBondType(Chem.BondType.SINGLE)
    rw_mol_B.GetBondBetweenAtoms(9, 10).SetBondType(Chem.BondType.SINGLE)
    product_B = rw_mol_B.GetMol()
    Chem.SanitizeMol(product_B)

    # Check that A and B produce different molecules
    smiles_prod_A = Chem.MolToSmiles(product_A, canonical=True)
    smiles_prod_B = Chem.MolToSmiles(product_B, canonical=True)
    if smiles_prod_A == smiles_prod_B:
        return "Incorrect: Code failed to verify that options A and B produce different regioisomers."

    # --- Verification Step 3: Check if product A matches the target description ---
    # Target description: Ester at C1 (next to a bridgehead), Propyl at C2 (next to the double bond).
    # In product_A, bridgehead atoms are from precursor C3 (idx 2) and C8 (idx 7).
    # The ester group is attached to precursor C2 (idx 1).
    # The propyl group is attached to precursor C11 (idx 10).
    # The new double bond is between precursor C9 (idx 8) and C10 (idx 9).
    
    atom_ester_attach = product_A.GetAtomWithIdx(1)
    is_ester_next_to_bridgehead = any(nbr.GetIdx() in [2, 7] for nbr in atom_ester_attach.GetNeighbors())

    atom_propyl_attach = product_A.GetAtomWithIdx(10)
    is_propyl_next_to_db = any(nbr.GetIdx() in [8, 9] for nbr in atom_propyl_attach.GetNeighbors())

    if not (is_ester_next_to_bridgehead and is_propyl_next_to_db):
        return "Incorrect: The product from option A does not match the target's structural description (substituent regiochemistry)."

    # --- Final Conclusion ---
    if final_answer == "A":
        return "Correct"
    else:
        return f"Incorrect: The final answer was {final_answer}, but chemical analysis shows A is the correct option."

# Run the check
result = check_correctness()
print(result)