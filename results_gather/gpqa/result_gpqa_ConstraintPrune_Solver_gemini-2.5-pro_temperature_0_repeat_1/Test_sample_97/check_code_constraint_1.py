# You may need to install the RDKit library: pip install rdkit-pypi
import sys
try:
    from rdkit import Chem
except ImportError:
    # This message will be shown if the user runs the code without rdkit
    print("Error: RDKit is not installed. Please install it by running 'pip install rdkit-pypi' in your terminal.", file=sys.stderr)
    sys.exit(1)

def check_reaction_correctness():
    """
    Checks the correctness of the identified starting material for the given reaction.
    Reaction: A + 1-propene --(Ru catalyst)--> 1-(prop-1-en-1-yl)-2-vinylcyclopentane
    """
    llm_provided_answer = "B"

    candidates = {
        "A": {"name": "2-methylbicyclo[3.1.0]hex-2-ene", "smiles": "CC1=CC2C1C2"},
        "B": {"name": "bicyclo[3.2.0]hept-6-ene", "smiles": "C1=CC2CCC1C2"},
        "C": {"name": "1,2-dimethylenecyclopentane", "smiles": "C=C1CCCC1=C"},
        "D": {"name": "2-methyl-3-methylenebicyclo[2.1.0]pentane", "smiles": "CC1C(=C)C2C12"}
    }

    product_smiles = "C/C=C/C1CCCC1C=C"  # 1-(prop-1-en-1-yl)-2-vinylcyclopentane
    coreactant_smiles = "C=CC"  # 1-propene

    # --- Constraint 1: Carbon Count Conservation ---
    product_mol = Chem.MolFromSmiles(product_smiles)
    coreactant_mol = Chem.MolFromSmiles(coreactant_smiles)
    
    # Use a query for carbon atoms for precision
    carbon_smart = Chem.MolFromSmarts("[#6]")
    product_carbons = len(product_mol.GetSubstructMatches(carbon_smart))
    coreactant_carbons = len(coreactant_mol.GetSubstructMatches(carbon_smart))
    
    required_carbons_in_A = product_carbons - coreactant_carbons
    
    if required_carbons_in_A != 7:
        return f"Internal logic error: Carbon count calculation is incorrect. Product has {product_carbons}C, co-reactant has {coreactant_carbons}C, requiring {required_carbons_in_A}C for A."

    passing_c1 = []
    for key, data in candidates.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        candidate_carbons = len(mol.GetSubstructMatches(carbon_smart))
        if candidate_carbons == required_carbons_in_A:
            passing_c1.append(key)
    
    # Note: All four candidates have 7 carbons, so this constraint does not prune any.
    # The LLM's reasoning that candidate A has 6 carbons was incorrect.

    # --- Constraint 2: Must be Bicyclic for ROCM ---
    passing_c2 = []
    for key in passing_c1:
        mol = Chem.MolFromSmiles(candidates[key]["smiles"])
        # GetRingInfo().NumRings() correctly identifies bicyclic systems.
        if mol.GetRingInfo().NumRings() > 1:
            passing_c2.append(key)
    
    # This step correctly prunes 'C' (1,2-dimethylenecyclopentane), which is monocyclic.
    if "C" in passing_c2:
        return "Reasoning Error: Candidate C is monocyclic and should have been pruned."
    
    # --- Constraint 3a: Must Contain a 5-Membered Ring to form the product core ---
    passing_c3a = []
    for key in passing_c2: # Check remaining candidates: A, B, D
        mol = Chem.MolFromSmiles(candidates[key]["smiles"])
        has_5_membered_ring = any(len(ring) == 5 for ring in mol.GetRingInfo().AtomRings())
        if has_5_membered_ring:
            passing_c3a.append(key)

    # This step correctly prunes 'D' (bicyclo[2.1.0]pentane derivative), which has no 5-membered ring.
    if "D" in passing_c3a:
        return "Reasoning Error: Candidate D does not contain a 5-membered ring and should have been pruned."

    # --- Constraint 3b: Reactive Double Bond Must NOT be in the 5-Membered Ring ---
    final_candidates = []
    for key in passing_c3a: # Check remaining candidates: A, B
        mol = Chem.MolFromSmiles(candidates[key]["smiles"])
        five_membered_rings_atoms = [set(ring) for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5]
        
        # Find the endocyclic double bond, which is the site of ROCM
        is_valid_candidate = True
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.IsInRing():
                atom1_idx, atom2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                # Check if this double bond is part of the 5-membered ring
                for ring_atoms in five_membered_rings_atoms:
                    if atom1_idx in ring_atoms and atom2_idx in ring_atoms:
                        # The double bond is in the 5-membered ring. ROCM would open this ring.
                        is_valid_candidate = False
                        break
                if not is_valid_candidate:
                    break
        
        if is_valid_candidate:
            final_candidates.append(key)

    # This step correctly prunes 'A', as its double bond is in the 5-membered ring.
    # It keeps 'B', where the double bond is in the 4-membered ring.
    if "A" in final_candidates:
        return "Reasoning Error: Candidate A should be pruned because its double bond is in the 5-membered ring, which would not be preserved."
    if "B" not in final_candidates:
        return "Reasoning Error: Candidate B was incorrectly pruned. It is the correct substrate."

    # --- Final Verification ---
    if len(final_candidates) == 1 and final_candidates[0] == llm_provided_answer:
        return "Correct"
    elif len(final_candidates) == 1:
        return f"Incorrect. The correct answer is {final_candidates[0]}, but the provided answer was {llm_provided_answer}."
    else:
        return f"Failed to determine a single correct answer. Remaining candidates: {final_candidates}"

# Run the check and print the result
try:
    result = check_reaction_correctness()
    print(result)
except SystemExit as e:
    # This handles the case where rdkit is not installed, preventing a crash.
    pass
except Exception as e:
    print(f"An unexpected error occurred: {e}", file=sys.stderr)
