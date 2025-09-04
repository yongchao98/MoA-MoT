# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def get_mol_from_name(name):
    """Converts an IUPAC name to an RDKit molecule object via a SMILES dictionary."""
    # RDKit does not have a built-in IUPAC name parser.
    # We use a pre-defined dictionary of SMILES strings for the compounds in the problem.
    name_to_smiles = {
        # Reactants and Products from hints and reactions
        "1,2-dimethyl-4-(propan-2-ylidene)cyclopentane": "CC1C(C)CC(=C(C)C)C1",
        "3,4-dimethylcyclopentan-1-one": "CC1CC(=O)C(C)C1",
        "1-cyano-3,4-dimethylcyclopentan-1-ol": "CC1CC(O)(C#N)C(C)C1",
        "1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol": "CC1CC(O)(CN)C(C)C1",
        # Options
        "2,2,3,4-tetramethylcyclobutan-1-one": "CC1(C)C(C)C(C)C1=O",
        "4-methylcycloheptan-1-one": "CC1CCCC(=O)CC1",
        "3,4-dimethylcyclohexan-1-one": "CC1CCC(=O)C(C)C1",
        "2,3,4-trimethylcyclopentan-1-one": "CC1C(C)C(C)C(=O)C1"
    }
    smiles = name_to_smiles.get(name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.AddHs(mol) if mol else None
    return None

def get_ketone_ring_size(mol):
    """Finds the ring size of a cyclic ketone."""
    if not mol: return 0
    # Find the carbon atom of the carbonyl group
    carbonyl_match = mol.GetSubstructMatch(Chem.MolFromSmarts('C=O'))
    if not carbonyl_match: return 0
    carbonyl_carbon_idx = carbonyl_match[0]

    # Check which ring contains this carbon
    for ring in Chem.GetSymmSSSR(mol):
        if carbonyl_carbon_idx in ring:
            return len(ring)
    return 0 # Not a cyclic ketone

def check_correctness():
    """
    Checks the correctness of the proposed answer by simulating the reaction scheme.
    """
    # The proposed answer from the LLM is C) 3,4-dimethylcyclohexan-1-one
    proposed_E_name = "3,4-dimethylcyclohexan-1-one"
    proposed_E_mol = get_mol_from_name(proposed_E_name)

    # --- Step 1: Identify Compound A from Hint (a) ---
    # The hint describes a Wittig reaction. We can identify A by a retro-Wittig on the product.
    wittig_product = get_mol_from_name("1,2-dimethyl-4-(propan-2-ylidene)cyclopentane")
    # The retro-reaction (ozonolysis) splits the C=C bond into two C=O groups.
    rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=[O].[C:2]=[O]')
    reactants = (wittig_product,)
    product_sets = rxn.RunReactants(reactants)
    
    ketones = [Chem.MolToSmiles(p, canonical=True) for p in product_sets[0]]
    expected_ketones = {
        Chem.MolToSmiles(get_mol_from_name("3,4-dimethylcyclopentan-1-one"), canonical=True),
        Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)C"), canonical=True) # Acetone
    }
    if set(ketones) != expected_ketones:
        return "Incorrect: Hint (a) implies that the Wittig reaction should start from 3,4-dimethylcyclopentan-1-one and acetone. The retro-synthesis did not yield these compounds."
    
    mol_A = get_mol_from_name("3,4-dimethylcyclopentan-1-one")

    # --- Step 2: Verify IR spectrum of A from Hint (b) ---
    ring_size_A = get_ketone_ring_size(mol_A)
    # 1750 cm-1 is characteristic of a strained 5-membered ring ketone.
    if ring_size_A != 5:
        return f"Incorrect: Compound A is expected to be a 5-membered ring ketone based on its IR spectrum (~1750 cm-1), but the derived structure has a {ring_size_A}-membered ring."

    # --- Step 3: A -> B (Cyanohydrin formation) ---
    rxn_A_to_B = AllChem.ReactionFromSmarts('[C:1]=[O:2]>>[C:1]([O:2])(C#N)')
    mol_B = rxn_A_to_B.RunReactants((mol_A,))[0][0]
    Chem.SanitizeMol(mol_B)
    expected_B = get_mol_from_name("1-cyano-3,4-dimethylcyclopentan-1-ol")
    if Chem.MolToSmiles(mol_B, canonical=True) != Chem.MolToSmiles(expected_B, canonical=True):
        return "Incorrect: The reaction of Compound A with HCN did not produce the expected cyanohydrin (Compound B)."

    # --- Step 4: B -> C (Nitrile reduction) ---
    rxn_B_to_C = AllChem.ReactionFromSmarts('[C:1]#[N:2]>>[CH2:1][NH2:2]')
    mol_C = rxn_B_to_C.RunReactants((mol_B,))[0][0]
    Chem.SanitizeMol(mol_C)
    expected_C = get_mol_from_name("1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol")
    if Chem.MolToSmiles(mol_C, canonical=True) != Chem.MolToSmiles(expected_C, canonical=True):
        return "Incorrect: The reduction of Compound B did not produce the expected amino alcohol (Compound C)."

    # --- Step 5: C -> E (Tiffeneau-Demjanov rearrangement) ---
    # This reaction involves ring expansion. The 5-membered ring in C should become a 6-membered ring in E.
    # The overall reaction is: C + HNO2 -> E + N2 + 2*H2O
    # We check atom conservation and ring size change.
    formula_C = Descriptors.rdMolDescriptors.CalcMolFormula(mol_C)
    formula_E = Descriptors.rdMolDescriptors.CalcMolFormula(proposed_E_mol)
    
    # Expected formulas based on the reaction C8H17NO -> C8H14O
    if formula_C != "C8H17NO":
        return f"Incorrect: The derived formula for Compound C is {formula_C}, but it should be C8H17NO."
    if formula_E != "C8H14O":
        return f"Incorrect: The proposed Compound E ({proposed_E_name}) has a formula of {formula_E}, but the product of the rearrangement should be C8H14O."

    ring_size_C = get_ketone_ring_size(mol_A) # C is not a ketone, but has the same ring as A
    ring_size_E = get_ketone_ring_size(proposed_E_mol)
    if not (ring_size_C == 5 and ring_size_E == 6):
        return f"Incorrect: The Tiffeneau-Demjanov rearrangement requires a ring expansion from 5 to 6 members. The check shows a transition from a {ring_size_C}-membered ring to a {ring_size_E}-membered ring, which is inconsistent with the proposed answer."

    # --- Step 6: Verify IR spectrum of E from Hint (b) ---
    # 1715 cm-1 is characteristic of a less-strained 6-membered ring ketone.
    if ring_size_E != 6:
        return f"Incorrect: The proposed Compound E is expected to be a 6-membered ring ketone based on its IR spectrum (~1715 cm-1), but it has a {ring_size_E}-membered ring."

    # --- All checks passed ---
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)