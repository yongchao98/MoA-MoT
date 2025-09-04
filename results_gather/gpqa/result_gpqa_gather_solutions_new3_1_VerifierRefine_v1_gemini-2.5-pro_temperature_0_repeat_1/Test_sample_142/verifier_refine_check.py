from rdkit import Chem
from rdkit.Chem import Descriptors

def get_formula(mol):
    """Calculates the molecular formula from an RDKit Mol object."""
    if not mol:
        return None
    return Descriptors.rdMolDescriptors.CalcMolFormula(mol)

def get_expected_start_formula(product_mol):
    """Calculates the expected formula for a starting material in a dehydration reaction."""
    if not product_mol:
        return None
    # Create a "product + H2O" molecule for formula comparison
    product_plus_water_mol = Chem.MolFromSmiles(Chem.MolToSmiles(product_mol) + ".O")
    return get_formula(product_plus_water_mol)

def get_aliphatic_ring_sizes(mol):
    """Returns a sorted list of non-aromatic ring sizes in the molecule."""
    if not mol:
        return []
    aliphatic_rings = []
    for ring_atoms in mol.GetRingInfo().AtomRings():
        is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms)
        if not is_aromatic:
            aliphatic_rings.append(len(ring_atoms))
    return sorted(aliphatic_rings)

def run_check():
    """
    Checks the correctness of the answer by verifying atom conservation and key mechanistic features
    for the two Pinacol-Pinacolone rearrangement reactions.
    """
    # --- Data Definition ---
    compounds = {
        "smiles_prod_1": "O=C1CCCCC1(c1ccc(C)cc1)c1ccc(C)cc1",
        "smiles_A_cyclopent": "OC1(C(O)(c2ccc(C)cc2)c2ccc(C)cc2)CCCC1",
        "smiles_A_cyclohex": "OC1(C(O)(c2ccc(C)cc2)c2ccc(C)cc2)CCCCC1",
        "smiles_start_2": "CC(O)C(O)(c1ccc(C)cc1)C(=O)OC",
        "smiles_B_correct": "CC(=O)C(c1ccc(C)cc1)C(=O)OC",
        "smiles_B_incorrect": "O=CC(C)(c1ccc(C)cc1)C(=O)OC",
    }

    options = {
        "A": {"A": "smiles_A_cyclohex", "B": "smiles_B_incorrect"},
        "B": {"A": "smiles_A_cyclopent", "B": "smiles_B_correct"},
        "C": {"A": "smiles_A_cyclopent", "B": "smiles_B_incorrect"},
        "D": {"A": "smiles_A_cyclohex", "B": "smiles_B_correct"},
    }

    llm_answer = "B"
    selected_option_keys = options.get(llm_answer)

    # --- Check Reaction 1 ---
    mol_prod_1 = Chem.MolFromSmiles(compounds["smiles_prod_1"])
    mol_A = Chem.MolFromSmiles(compounds[selected_option_keys["A"]])
    
    expected_A_formula = get_expected_start_formula(mol_prod_1)
    actual_A_formula = get_formula(mol_A)
    
    if expected_A_formula != actual_A_formula:
        return (f"Incorrectness in Reaction 1: The starting material A from option {llm_answer} does not have the correct molecular formula for a dehydration reaction. "
                f"Product formula is {get_formula(mol_prod_1)}. Expected formula for A is {expected_A_formula}. "
                f"Actual formula for A is {actual_A_formula}.")

    start_A_rings = get_aliphatic_ring_sizes(mol_A)
    prod_1_rings = get_aliphatic_ring_sizes(mol_prod_1)

    if not (5 in start_A_rings and 6 in prod_1_rings):
        return (f"Incorrectness in Reaction 1: The mechanism requires ring expansion from a 5-membered ring to a 6-membered ring. "
                f"The chosen starting material A (aliphatic rings: {start_A_rings}) or product (aliphatic rings: {prod_1_rings}) does not fit this mechanism.")

    # --- Check Reaction 2 ---
    mol_start_2 = Chem.MolFromSmiles(compounds["smiles_start_2"])
    mol_B = Chem.MolFromSmiles(compounds[selected_option_keys["B"]])

    expected_start_2_formula = get_expected_start_formula(mol_B)
    actual_start_2_formula = get_formula(mol_start_2)

    if expected_start_2_formula != actual_start_2_formula:
        return (f"Incorrectness in Reaction 2: The product B from option {llm_answer} does not have the correct molecular formula for a dehydration reaction. "
                f"Starting material formula is {actual_start_2_formula}. Expected formula for start is {expected_start_2_formula}, derived from product B with formula {get_formula(mol_B)}.")

    return "Correct"

# Execute the check
result = run_check()
print(result)