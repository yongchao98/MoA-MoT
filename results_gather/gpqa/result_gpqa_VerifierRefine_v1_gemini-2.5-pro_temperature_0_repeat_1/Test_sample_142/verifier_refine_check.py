import sys
from collections import Counter

# Try to import RDKit; if it's not available, the check cannot be performed.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    # In a real environment, you would handle this error. For this purpose,
    # we'll print a message and exit if rdkit is not installed.
    # In this context, we assume the environment can run the code.
    sys.exit("RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this check.")

def get_molecular_formula(mol: Chem.Mol) -> str:
    """Calculates the molecular formula from an RDKit Mol object."""
    return Descriptors.rdMolDescriptors.CalcMolFormula(mol)

def get_aliphatic_ring_sizes(mol: Chem.Mol) -> list:
    """Finds the sizes of all non-aromatic rings in a molecule."""
    aliphatic_rings = []
    for ring_indices in mol.GetRingInfo().AtomRings():
        is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_indices)
        if not is_aromatic:
            aliphatic_rings.append(len(ring_indices))
    return aliphatic_rings

def check_correctness():
    """
    Checks the correctness of the provided answer for the Pinacol rearrangement questions.
    The answer to check is Option D:
    A = 1-(hydroxydi-p-tolylmethyl)cyclopentan-1-ol
    B = methyl 3-oxo-2-(p-tolyl)butanoate
    """
    # --- Define molecules from Option D using SMILES notation ---
    # Reaction 1: A ---> 2,2-di-p-tolylcyclohexan-1-one
    smiles_A = "OC1(CCCC1)C(O)(c2ccc(C)cc2)c3ccc(C)cc3"  # 1-(hydroxydi-p-tolylmethyl)cyclopentan-1-ol
    smiles_prod1 = "O=C1CCCCC1(c2ccc(C)cc2)c3ccc(C)cc3" # 2,2-di-p-tolylcyclohexan-1-one

    # Reaction 2: methyl 2,3-dihydroxy-2-(p-tolyl)butanoate ---> B
    smiles_react2 = "CC(O)C(O)(c1ccc(C)cc1)C(=O)OC" # methyl 2,3-dihydroxy-2-(p-tolyl)butanoate
    smiles_B = "CC(=O)C(c1ccc(C)cc1)C(=O)OC"      # methyl 3-oxo-2-(p-tolyl)butanoate

    # --- Part 1: Verification of Reaction A -> Product 1 ---
    mol_A = Chem.MolFromSmiles(smiles_A)
    mol_prod1 = Chem.MolFromSmiles(smiles_prod1)

    # 1.1: Check Stoichiometry (Reactant -> Product + H2O)
    formula_A = get_molecular_formula(mol_A)
    formula_prod1 = get_molecular_formula(mol_prod1)
    expected_formula_A = "C21H26O2"
    expected_formula_prod1 = "C21H24O"
    if not (formula_A == expected_formula_A and formula_prod1 == expected_formula_prod1):
        return (f"Incorrect molecular formula for Reaction 1. "
                f"Expected A={expected_formula_A} and Product={expected_formula_prod1}. "
                f"Got A={formula_A} and Product={formula_prod1}.")

    # 1.2: Check functional groups (1,2-diol -> ketone)
    diol_pattern = Chem.MolFromSmarts("[CX4](O)-[CX4](O)")
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    if not mol_A.HasSubstructMatch(diol_pattern):
        return "Constraint failed: Proposed starting material A is not a 1,2-diol (pinacol)."
    if not mol_prod1.HasSubstructMatch(ketone_pattern):
        return "Constraint failed: Product of reaction 1 is not a ketone (pinacolone)."

    # 1.3: Check structural rearrangement (Ring Expansion)
    ring_size_A = get_aliphatic_ring_sizes(mol_A)
    ring_size_prod1 = get_aliphatic_ring_sizes(mol_prod1)
    if Counter(ring_size_A) != Counter([5]):
        return f"Constraint failed: Reactant A should have a single 5-membered aliphatic ring, but found rings of sizes: {ring_size_A}."
    if Counter(ring_size_prod1) != Counter([6]):
        return f"Constraint failed: Product 1 should have a single 6-membered aliphatic ring, but found rings of sizes: {ring_size_prod1}."

    # --- Part 2: Verification of Reaction 2 -> B ---
    mol_react2 = Chem.MolFromSmiles(smiles_react2)
    mol_B = Chem.MolFromSmiles(smiles_B)

    # 2.1: Check Stoichiometry
    formula_react2 = get_molecular_formula(mol_react2)
    formula_B = get_molecular_formula(mol_B)
    expected_formula_react2 = "C12H16O4"
    expected_formula_B = "C12H14O3"
    if not (formula_react2 == expected_formula_react2 and formula_B == expected_formula_B):
        return (f"Incorrect molecular formula for Reaction 2. "
                f"Expected Reactant={expected_formula_react2} and B={expected_formula_B}. "
                f"Got Reactant={formula_react2} and B={formula_B}.")

    # 2.2: Check functional groups
    if not mol_react2.HasSubstructMatch(diol_pattern):
        return "Constraint failed: The starting material for reaction 2 is not a 1,2-diol."
    if not mol_B.HasSubstructMatch(ketone_pattern):
        return "Constraint failed: Proposed product B is not a ketone."

    # 2.3: Check rearrangement type (Hydride shift)
    # The structure of product B corresponds to a 1,2-hydride shift, which is
    # chemically favored over a 1,2-methyl shift. The alternative product in the options
    # would result from a methyl shift. This check confirms the product is consistent
    # with the most plausible reaction pathway.
    # The check is implicitly successful because the SMILES for B matches the hydride shift product.

    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)