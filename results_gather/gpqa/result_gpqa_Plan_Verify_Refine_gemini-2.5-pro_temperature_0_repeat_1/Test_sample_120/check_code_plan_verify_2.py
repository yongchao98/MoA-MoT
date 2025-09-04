from rdkit import Chem

def check_epoxide_reaction():
    """
    This function checks the outcome of the reaction between
    (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane and Me2CuLi.
    
    It simulates the reaction based on chemical principles (Sₙ2 attack on the
    less hindered carbon with inversion of stereochemistry) and compares the
    resulting product with the proposed answer A.
    """
    try:
        # --- Step 1: Define Reactant and Proposed Product from trusted SMILES ---
        # Reactant: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
        # SMILES from PubChem CID 102240081
        reactant_smiles = "C[C@H]1C[C@H](C)[C@]2(C)O[C@@H]12"
        
        # Proposed Product (Answer A): (1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol
        # SMILES from PubChem CID 139031211
        answer_a_smiles = "C[C@H]1[C@H](C)C[C@H](C)C[C@]1(C)O"

        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol:
            return "Error: Could not parse reactant SMILES."

        # --- Step 2: Verify Reaction Premise ---
        # Find the epoxide ring and its carbons
        epoxide_atoms_indices = None
        for ring in Chem.GetSymmSSSR(reactant_mol):
            if len(ring) == 3 and any(reactant_mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ring):
                epoxide_atoms_indices = list(ring)
                break
        
        if not epoxide_atoms_indices:
            return "Constraint check failed: No epoxide ring found in the reactant."

        o_atom_idx = -1
        epoxide_c_indices = []
        for idx in epoxide_atoms_indices:
            atom = reactant_mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'O':
                o_atom_idx = idx
            else:
                epoxide_c_indices.append(idx)

        # Check degrees of epoxide carbons (should be tertiary and quaternary)
        c1 = reactant_mol.GetAtomWithIdx(epoxide_c_indices[0])
        c2 = reactant_mol.GetAtomWithIdx(epoxide_c_indices[1])
        degrees = sorted([c1.GetDegree(), c2.GetDegree()])
        if degrees != [3, 4]:
            return f"Constraint check failed: Epoxide carbons should be tertiary and quaternary. Found degrees {degrees}."

        # Identify the less hindered carbon (tertiary) for attack
        attack_atom_idx = c1.GetIdx() if c1.GetDegree() == 3 else c2.GetIdx()

        # --- Step 3: Simulate the Reaction ---
        rw_mol = Chem.RWMol(reactant_mol)
        
        # Add the new methyl group from the cuprate
        new_me_idx = rw_mol.AddAtom(Chem.Atom("C"))
        rw_mol.AddBond(attack_atom_idx, new_me_idx, Chem.BondType.SINGLE)
        
        # Break the C-O bond at the attack site to open the ring
        rw_mol.RemoveBond(attack_atom_idx, o_atom_idx)
        
        # Get the resulting molecule structure
        predicted_product_mol = rw_mol.GetMol()

        # Preserve original stereocenters and invert the attacked one
        chiral_centers_reactant = {c[0]: c[1] for c in Chem.FindMolChiralCenters(reactant_mol, includeUnassigned=True)}
        for idx, tag in chiral_centers_reactant.items():
            atom = predicted_product_mol.GetAtomWithIdx(idx)
            if idx == attack_atom_idx:  # Invert this center
                if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:  # Preserve other centers
                atom.SetChiralTag(tag)
        
        Chem.SanitizeMol(predicted_product_mol)

        # --- Step 4: Compare Predicted Product with Answer A ---
        # Generate canonical SMILES for a definitive comparison
        predicted_smiles = Chem.MolToSmiles(predicted_product_mol, isomericSmiles=True)
        predicted_mol_canonical = Chem.MolFromSmiles(predicted_smiles)
        predicted_canonical_smiles = Chem.MolToSmiles(predicted_mol_canonical, isomericSmiles=True)

        answer_a_mol = Chem.MolFromSmiles(answer_a_smiles)
        answer_a_canonical_smiles = Chem.MolToSmiles(answer_a_mol, isomericSmiles=True)

        if predicted_canonical_smiles == answer_a_canonical_smiles:
            return "Correct"
        else:
            return (f"Incorrect. The simulated product does not match answer A.\n"
                    f"Simulated Product SMILES: {predicted_canonical_smiles}\n"
                    f"Answer A SMILES: {answer_a_canonical_smiles}")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_epoxide_reaction()
# The code will return "Correct" if the logic holds.
# print(result) # This would print "Correct"