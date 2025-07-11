# To run this code, you need to have the RDKit library installed.
# You can install it using pip: pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def solve_chemistry_problem():
    """
    Solves the multi-step synthesis problem to identify Compound 3.
    """
    # Starting material: Terpinolene
    # IUPAC: 1-methyl-4-(propan-2-ylidene)cyclohex-1-ene
    terpinolene_smiles = 'CC1=CCC(C(=C(C)C))CC1'
    mol = Chem.MolFromSmiles(terpinolene_smiles)
    if not mol:
        print("Error: Could not parse Terpinolene SMILES.")
        return

    print(f"Starting Material: Terpinolene")
    print(f"SMILES: {Chem.MolToSmiles(mol)}\n")

    # --- Step 1: Epoxidation with m-CPBA ---
    # m-CPBA epoxidizes the more substituted double bond.
    # We need to find the most substituted C=C bond in terpinolene.
    # The exocyclic double bond is tetrasubstituted, while the endocyclic one is trisubstituted.
    # Therefore, the reaction occurs at the exocyclic double bond.
    
    # Let's find the exocyclic double bond atoms C=C(C)C
    # The SMARTS pattern for this specific substructure in terpinolene is '[#6](=C(C)C)-C'
    patt = Chem.MolFromSmarts('[#6;R]=[#6;R]') # endocyclic
    endo_bond_indices = mol.GetSubstructMatch(patt)
    
    patt = Chem.MolFromSmarts('[#6;R]-C(=C(-C)-C)') # exocyclic C=C(C)C attached to a ring carbon
    exo_bond_indices = mol.GetSubstructMatch(patt)
    
    # The reaction will happen on the atoms from the exocyclic pattern
    # The C=C atoms are at indices 1 and 2 of this pattern match
    c1_idx, c2_idx = exo_bond_indices[1], exo_bond_indices[2]

    # Build Compound 1 (epoxide) by manipulating the molecule
    emol = Chem.EditableMol(mol)
    emol.RemoveBond(c1_idx, c2_idx)
    o_idx = emol.AddAtom(Chem.Atom(8)) # Oxygen atom
    emol.AddBond(c1_idx, o_idx, Chem.BondType.SINGLE)
    emol.AddBond(c2_idx, o_idx, Chem.BondType.SINGLE)
    
    compound1 = emol.GetMol()
    Chem.SanitizeMol(compound1)
    
    print("--- Step 1: Epoxidation ---")
    print("Terpinolene reacts with m-CPBA. The more substituted exocyclic double bond is epoxidized.")
    print(f"Compound 1 (Epoxide) SMILES: {Chem.MolToSmiles(compound1)}\n")

    # --- Step 2: Formation of Thiirane ---
    # The epoxide is converted to a thiirane with N,N-dimethyl thioformamide and acid.
    # This means the Oxygen atom of the epoxide ring is replaced by a Sulfur atom.
    
    # Find the oxygen atom in the 3-membered ring and replace it with sulfur
    compound2 = Chem.RWMol(compound1)
    for atom in compound2.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.IsInRingSize(3):
            atom.SetAtomicNum(16) # 16 is Sulfur
            break
            
    Chem.SanitizeMol(compound2)
    
    print("--- Step 2: Thiirane Formation ---")
    print("Compound 1 reacts to replace the epoxide oxygen with sulfur.")
    print(f"Compound 2 (Thiirane) SMILES: {Chem.MolToSmiles(compound2)}\n")

    # --- Step 3: Reduction with LiAlH4 ---
    # LiAlH4 reduces the thiirane back to an alkene (desulfurization).
    # We need to find the sulfur atom, remove it, and form a double bond between its neighbors.
    
    emol3 = Chem.EditableMol(compound2)
    s_idx = -1
    neighbors = []
    for atom in compound2.GetAtoms():
        if atom.GetAtomicNum() == 16 and atom.IsInRingSize(3):
            s_idx = atom.GetIdx()
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            break
    
    # Remove the sulfur atom (this also removes its bonds)
    emol3.RemoveAtom(s_idx)
    # Add a double bond between the two carbons that were attached to the sulfur
    emol3.AddBond(neighbors[0], neighbors[1], Chem.BondType.DOUBLE)
    
    compound3 = emol3.GetMol()
    Chem.SanitizeMol(compound3)

    print("--- Step 3: Reduction ---")
    print("Compound 2 is reduced with LiAlH4, which removes the sulfur and reforms the double bond.")
    print(f"Compound 3 SMILES: {Chem.MolToSmiles(compound3)}\n")

    # --- Final Answer ---
    print("--- Final Product Identification ---")
    # As we can see, the final product is the same as the starting material.
    final_mol = compound3
    formula = Descriptors.rdMolDescriptors.CalcMolFormula(final_mol)
    
    # Count atoms for the "output each number" request
    num_c = 0
    num_h = 0
    for atom in final_mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            num_c += 1
        # Calculate hydrogens, which are implicit
    # A more direct way to get H count
    num_h = Descriptors.rdMolDescriptors.CalcNumHs(final_mol)

    print("Compound 3 is the starting material, Terpinolene.")
    print(f"Molecular Formula: {formula}")
    print(f"The numbers in the final formula are: Carbon = {num_c}, Hydrogen = {num_h}")


if __name__ == '__main__':
    solve_chemistry_problem()
