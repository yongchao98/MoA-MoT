import sys
try:
    from rdkit import Chem
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def check_michael_addition_answer():
    """
    Checks the correctness of the given answer for three Michael addition reactions.
    It verifies the structures of products A, B, and reactant C.
    """
    # --- Part 1: Analysis of Reaction A ---
    # Reaction: dimethyl malonate + methyl (E)-3-(p-tolyl)acrylate -> A
    # Michael Donor: Enolate of dimethyl malonate [-CH(COOMe)2]
    # Michael Acceptor: methyl (E)-3-(p-tolyl)acrylate [p-tolyl-CH(beta)=CH(alpha)-COOMe]
    # The malonate enolate attacks the beta-carbon. The resulting structure is:
    # p-tolyl-CH(CH(COOMe)2)-CH2-COOMe
    # This corresponds to the IUPAC name: trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate.
    # We generate a canonical SMILES string for this expected structure.
    expected_A_smiles = "COC(=O)CC(c1ccc(C)cc1)C(C(=O)OC)C(=O)OC"
    mol_expected_A = Chem.MolFromSmiles(expected_A_smiles)
    canonical_expected_A = Chem.MolToSmiles(mol_expected_A, canonical=True)

    # --- Part 2: Analysis of Reaction B ---
    # Reaction: 1-(cyclohex-1-en-1-yl)piperidine + (E)-but-2-enenitrile -> B
    # This is a Stork enamine alkylation, a variant of the Michael reaction.
    # Michael Donor: Enamine of cyclohexanone. The nucleophilic carbon is C-alpha to the original carbonyl.
    # Michael Acceptor: (E)-but-2-enenitrile [CH3-CH(beta)=CH(alpha)-CN]
    # The enamine's alpha-carbon attacks the beta-carbon of the nitrile.
    # Subsequent acidic workup (H3O+) hydrolyzes the intermediate iminium salt back to a ketone.
    # The final product is cyclohexanone with a -CH(CH3)CH2CN group at the alpha-position (C2).
    # This corresponds to the IUPAC name: 3-(2-oxocyclohexyl)butanenitrile.
    # The keto form is the major, stable product.
    expected_B_smiles = "CC(CC#N)C1CCCCC1=O"
    mol_expected_B = Chem.MolFromSmiles(expected_B_smiles)
    canonical_expected_B = Chem.MolToSmiles(mol_expected_B, canonical=True)

    # --- Part 3: Analysis of Reaction C (Retrosynthesis) ---
    # Reaction: C + but-3-en-2-one -> 2-(3-oxobutyl)cyclohexane-1,3-dione
    # Product: 2-(3-oxobutyl)cyclohexane-1,3-dione
    # Michael Acceptor: but-3-en-2-one (MVK). In the product, this becomes the "3-oxobutyl" group [-CH2-CH2-CO-CH3].
    # The bond is formed between the beta-carbon of MVK and the nucleophilic carbon of the donor C.
    # The "3-oxobutyl" group is attached to C2 of the cyclohexane-1,3-dione ring.
    # Therefore, C2 must be the nucleophilic carbon. This carbon is highly acidic as it's between two carbonyls.
    # The Michael donor (reactant C) must be cyclohexane-1,3-dione.
    expected_C_smiles = "O=C1CC(=O)CC1"
    mol_expected_C = Chem.MolFromSmiles(expected_C_smiles)
    canonical_expected_C = Chem.MolToSmiles(mol_expected_C, canonical=True)

    # --- Verification of the provided answer (Option C) ---
    # The answer claims:
    # A = trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate
    # B = 3-(2-oxocyclohexyl)butanenitrile
    # C = cyclohexane-1,3-dione
    
    # We convert these names to SMILES to compare with our expected structures.
    answer_A_smiles = "COC(=O)CC(c1ccc(C)cc1)C(C(=O)OC)C(=O)OC"
    answer_B_smiles = "CC(CC#N)C1CCCCC1=O"
    answer_C_smiles = "O=C1CC(=O)CC1"

    # Check A
    mol_answer_A = Chem.MolFromSmiles(answer_A_smiles)
    canonical_answer_A = Chem.MolToSmiles(mol_answer_A, canonical=True)
    if canonical_answer_A != canonical_expected_A:
        return (f"Incorrect product A. The provided answer corresponds to the structure {canonical_answer_A}, "
                f"but the correct product of the reaction is {canonical_expected_A}. "
                "The p-tolyl group should be on C2 of the propane chain, not C3.")

    # Check B
    mol_answer_B = Chem.MolFromSmiles(answer_B_smiles)
    canonical_answer_B = Chem.MolToSmiles(mol_answer_B, canonical=True)
    if canonical_answer_B != canonical_expected_B:
        return (f"Incorrect product B. The provided answer corresponds to the structure {canonical_answer_B}, "
                f"but the correct product is {canonical_expected_B}. "
                "The question asks for the major final product, which is the keto form (2-oxocyclohexyl), not an enol tautomer.")

    # Check C
    mol_answer_C = Chem.MolFromSmiles(answer_C_smiles)
    canonical_answer_C = Chem.MolToSmiles(mol_answer_C, canonical=True)
    if canonical_answer_C != canonical_expected_C:
        return (f"Incorrect reactant C. The provided answer corresponds to the structure {canonical_answer_C}, "
                f"but the correct Michael donor is {canonical_expected_C}. "
                "The reactant is the dione itself, which is then deprotonated in situ to form the nucleophile.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_michael_addition_answer()
print(result)