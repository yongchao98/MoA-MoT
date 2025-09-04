from rdkit import Chem

def get_stereochem_from_smiles(smiles_string):
    """
    Parses a SMILES string and returns the CIP stereochemical labels (R/S)
    for the relevant stereocenters in the given problem.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        return "Invalid SMILES"
    
    # Find the atoms corresponding to C1, C2, and C4 of the cyclohexane ring
    # This requires finding the substructure match
    pattern = Chem.MolFromSmarts('C1(C)(O)C(OC)CC(C(C)C)C1')
    if not mol.HasSubstructMatch(pattern):
        # Handle ester variant
        pattern = Chem.MolFromSmarts('C1(C)(OC(=O)CC)C(OC)CC(C(C)C)C1')
        if not mol.HasSubstructMatch(pattern):
            return "Structure does not match expected backbone"

    match_indices = mol.GetSubstructMatch(pattern)
    
    # Map pattern atoms to molecule atoms
    # Pattern: C1(idx 0), C2(idx 2), C4(idx 5)
    c1_idx, c2_idx, c4_idx = match_indices[0], match_indices[2], match_indices[5]

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    stereo_info = {}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in [c1_idx, c2_idx, c4_idx] and atom.HasProp('_CIPCode'):
            if atom.GetIdx() == c1_idx:
                stereo_info['C1'] = atom.GetProp('_CIPCode')
            elif atom.GetIdx() == c2_idx:
                stereo_info['C2'] = atom.GetProp('_CIPCode')
            elif atom.GetIdx() == c4_idx:
                stereo_info['C4'] = atom.GetProp('_CIPCode')
                
    return f"C1:{stereo_info.get('C1', 'N/A')}, C2:{stereo_info.get('C2', 'N/A')}, C4:{stereo_info.get('C4', 'N/A')}"


# --- Step-by-step analysis based on chemical principles ---

# Step 1: Hydrogenation of (R)-(+)-Limonene
# (R)-Limonene: C[C@@H]1CC=C(C(=C)C)CC1
# Pd/C selectively hydrogenates the less substituted exocyclic double bond.
# Product 1 is (R)-4-isopropyl-1-methylcyclohex-1-ene. The (R) stereocenter at C4 is unaffected.
# p1_smiles = 'C[C@@H]1CC=C(C)C(C(C)C)C1'

# Step 2: Epoxidation of Product 1
# mCPBA epoxidation is stereoselective. The bulky (R)-isopropyl group is on the 'wedge' face.
# The reagent attacks from the less hindered 'dash' face, creating an epoxide that is trans to the isopropyl group.
# This process creates two new stereocenters. The resulting epoxide is (1R, 2S, 4R)-1,2-epoxy-4-isopropyl-1-methylcyclohexane.
# p2_smiles = 'C[C@@]1(O2)[C@@H](C[C@H](C(C)C)CC1)C2'

# Step 3: Epoxide Ring-Opening with Sodium Methoxide (NaOMe)
# NaOMe is a strong nucleophile, and the reaction proceeds under basic conditions.
# The nucleophilic attack (SN2) occurs at the less sterically hindered carbon of the epoxide, which is C2 (secondary) rather than C1 (tertiary).
# SN2 attack proceeds with inversion of configuration at the attacked center. C2 was (S), so it becomes (R).
# The stereocenter that is not attacked, C1, retains its configuration. C1 was (R), so it remains (R).
# The C4 stereocenter is also unaffected and remains (R).
# Therefore, the resulting alcohol (Product 3) must be the (1R, 2R, 4R) stereoisomer.
predicted_p3_smiles = 'C[C@@]1(O)[C@@H](OC)CC[C@H](C(C)C)C1'

# Step 4: Esterification
# Steglich esterification converts the alcohol to an ester without affecting any stereocenters.
# The final product (Product 4) must therefore be (1R, 2R, 4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate.
predicted_p4_smiles = 'C[C@@]1(OC(=O)CC)[C@@H](OC)CC[C@H](C(C)C)C1'

# --- Verification ---

# The LLM's answer is A, which is named (1S,2R,4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate.
# Let's generate the SMILES for option A to compare.
option_A_smiles = 'C[C@]1(OC(=O)CC)[C@@H](OC)CC[C@H](C(C)C)C1'

# Get the canonical stereochemistry for both the predicted product and the LLM's answer.
predicted_stereo = get_stereochem_from_smiles(predicted_p4_smiles)
llm_answer_stereo = get_stereochem_from_smiles(option_A_smiles)

if predicted_stereo == llm_answer_stereo:
    print("Correct")
else:
    reason = f"""The provided answer is incorrect.

A step-by-step analysis based on established chemical principles leads to a different product than the one selected.

1.  **Predicted Product Stereochemistry**: The reaction sequence should yield the (1R, 2R, 4R) stereoisomer.
    *   The key step is the basic ring-opening of the (1R, 2S, 4R)-epoxide.
    *   The methoxide nucleophile attacks the less-hindered C2, inverting its configuration from (S) to (R).
    *   The configuration at C1 is not affected by this reaction and remains (R).
    *   The predicted stereochemistry is therefore (1R, 2R, 4R).
    *   My analysis results in: {predicted_stereo}.

2.  **LLM's Answer Stereochemistry**: The selected answer, option A, is the (1S, 2R, 4R) stereoisomer.
    *   The LLM's answer corresponds to: {llm_answer_stereo}.

3.  **Discrepancy**: The conflict is at the C1 stereocenter. The LLM's answer claims it is (S), while chemical principles dictate it should be (R). The LLM's written explanation correctly identifies the reaction pathway (attack at C2) but then incorrectly determines the resulting stereochemistry at C1, a classic reasoning error. The isomer in option A is not the expected product of this reaction.
"""
    print(reason)
