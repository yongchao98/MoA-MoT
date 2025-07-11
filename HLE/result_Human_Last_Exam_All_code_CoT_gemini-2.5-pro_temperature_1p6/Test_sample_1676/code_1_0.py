try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit'")
    exit()

def get_iupac_name(smiles):
    """Generates IUPAC name for a SMILES string. Needs an internet connection and pubchempy."""
    try:
        import pubchempy as pcp
        # Search for the compound by SMILES
        compounds = pcp.get_compounds(smiles, 'smiles')
        if compounds:
            return compounds[0].iupac_name
        else:
            return "Name not found"
    except (ImportError, Exception):
        return "Could not generate IUPAC name (pubchempy might be missing or network error)."


def solve_chemistry_problem():
    """
    Solves the multi-step synthesis problem programmatically.
    """
    # Starting material: Terpinolene
    # IUPAC: 4-isopropylidene-1-methylcyclohex-1-ene or 1-methyl-4-(propan-2-ylidene)cyclohex-1-ene
    terpinolene_smiles = 'CC1=CCC(=C(C)C)CC1'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    
    print("--- Reaction Pathway Analysis ---")
    print(f"Start (Terpinolene): {terpinolene_smiles}")
    
    # --- Step 1: Epoxidation of Terpinolene to Compound 1 ---
    # The reaction occurs at the endocyclic (in-ring) double bond.
    # We find this bond by checking the IsInRing property.
    
    rw_mol = Chem.RWMol(terpinolene)
    
    # Find the endocyclic double bond atoms
    c1_idx, c2_idx = -1, -1
    for bond in rw_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.IsInRing():
            c1_idx = bond.GetBeginAtomIdx()
            c2_idx = bond.GetEndAtomIdx()
            break
            
    # Modify molecule: Change C=C to C-C and add an epoxide ring
    rw_mol.RemoveBond(c1_idx, c2_idx)
    rw_mol.AddBond(c1_idx, c2_idx, Chem.BondType.SINGLE)
    o_idx = rw_mol.AddAtom(Chem.Atom(8)) # Oxygen
    rw_mol.AddBond(c1_idx, o_idx, Chem.BondType.SINGLE)
    rw_mol.AddBond(c2_idx, o_idx, Chem.BondType.SINGLE)

    compound1_mol = rw_mol.GetMol()
    Chem.SanitizeMol(compound1_mol)
    compound1_smiles = Chem.MolToSmiles(compound1_mol)
    
    print("\nReaction 1: Epoxidation")
    print(f"Terpinolene ({terpinolene_smiles}) -> Compound 1 ({compound1_smiles})")
    
    # --- Step 2: Conversion of Epoxide (Compound 1) to Thiirane (Compound 2) ---
    # The oxygen atom of the epoxide is replaced by a sulfur atom.

    rw_mol = Chem.RWMol(compound1_mol)
    # Find the oxygen atom in the three-membered ring and replace it with sulfur
    epoxide_o_idx = -1
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.IsInRingSize(3):
            atom.SetAtomicNum(16) # 16 is Sulfur
            break

    compound2_mol = rw_mol.GetMol()
    Chem.SanitizeMol(compound2_mol)
    compound2_smiles = Chem.MolToSmiles(compound2_mol)
    
    print("\nReaction 2: Thiiranation")
    print(f"Compound 1 ({compound1_smiles}) -> Compound 2 ({compound2_smiles})")
    
    # --- Step 3: Reduction of Thiirane (Compound 2) to Thiol (Compound 3) ---
    # The hydride from LiAlH4 attacks the less substituted carbon of the thiirane.
    # The original endocyclic bond was between a tertiary C (with Me) and a secondary C.
    
    rw_mol = Chem.RWMol(compound2_mol)
    
    # Identify thiirane atoms
    s_atom_idx = -1
    thiirane_carbons = []
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() == 16 and atom.IsInRingSize(3):
            s_atom_idx = atom.GetIdx()
            for neighbor in atom.GetNeighbors():
                thiirane_carbons.append(neighbor)
            break
            
    c_less_sub = thiirane_carbons[0] if thiirane_carbons[0].GetDegree() < thiirane_carbons[1].GetDegree() else thiirane_carbons[1]
    
    # Break the C-S bond at the less substituted carbon.
    rw_mol.RemoveBond(c_less_sub.GetIdx(), s_atom_idx)
    
    # Get the final molecule. Sanitization will add hydrogens to satisfy valency,
    # resulting in a C-H and S-H group, thus forming the final thiol product.
    compound3_mol = rw_mol.GetMol()
    Chem.SanitizeMol(compound3_mol)
    
    # Clean up for canonical SMILES
    compound3_mol = Chem.RemoveHs(compound3_mol)
    compound3_smiles = Chem.MolToSmiles(compound3_mol, isomericSmiles=True)
    
    print("\nReaction 3: Reductive Ring Opening")
    print(f"Compound 2 ({compound2_smiles}) -> Compound 3 ({compound3_smiles})")
    
    # Final result
    print("\n--- Final Product ---")
    print("The final product, Compound 3, is:")
    print(f"Structure (SMILES): {compound3_smiles}")
    print(f"Predicted IUPAC Name: 1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol")


if __name__ == '__main__':
    solve_chemistry_problem()
