import sys
from typing import Dict

# Install rdkit if not available
try:
    from rdkit import Chem
except ImportError:
    print("RDKit not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def get_diester_config_from_smiles(smiles: str) -> str:
    """
    Determines the diester configuration (cis/trans) from a SMILES string using RDKit.
    It finds the two chiral carbons attached to ester groups and compares their
    stereochemical tags. Same tags imply a 'cis' relationship for adjacent
    substituents in a ring, while different tags imply 'trans'.
    
    Args:
        smiles: The SMILES string of the molecule.
        
    Returns:
        'cis' or 'trans' string indicating the diester configuration.
        
    Raises:
        ValueError: If the SMILES is invalid or the expected structure is not found.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # Use a SMARTS pattern to find chiral carbons attached to an ester group.
    # Pattern: [Chiral Carbon with 1 Hydrogen](-Carbonyl-Oxygen-Carbon)
    patt = Chem.MolFromSmarts('[C;H1;@,@@](C(=O)OC)')
    matches = mol.GetSubstructMatches(patt)
    
    if len(matches) != 2:
        raise ValueError("Could not find exactly two chiral carbons with ester groups.")

    # Get the indices of the two chiral carbons.
    idx1, idx2 = matches[0][0], matches[1][0]
    
    # Verify that these carbons are bonded to each other.
    bond = mol.GetBondBetweenAtoms(idx1, idx2)
    if not bond:
        raise ValueError("Ester-bearing carbons are not adjacent.")

    # Get the chiral tags for each atom.
    # RDKit tags: 1=CHI_TETRAHEDRAL_CW ('@'), 2=CHI_TETRAHEDRAL_CCW ('@@')
    tag1 = mol.GetAtomWithIdx(idx1).GetChiralTag()
    tag2 = mol.GetAtomWithIdx(idx2).GetChiralTag()
    
    # If the tags are the same (e.g., both CW or both CCW), the relationship is cis.
    # If they are different, the relationship is trans.
    if tag1 == tag2:
        return "cis"
    else:
        return "trans"

def check_answer():
    """
    Checks the correctness of the provided LLM answer by verifying its premises.
    """
    # SMILES strings for each option from the question
    options_smiles = {
        'A': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O',
        'B': 'O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O',
        'C': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O',
        'D': 'O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O'
    }

    # Stereochemical claims made in the provided LLM answer's code
    llm_claims = {
        'A': {'diester_config': 'trans'},
        'B': {'diester_config': 'cis'},
        'C': {'diester_config': 'cis'},
        'D': {'diester_config': 'trans'}  # This is the claim to be scrutinized
    }

    errors = []
    for option, claimed_props in llm_claims.items():
        smiles = options_smiles[option]
        claimed_config = claimed_props['diester_config']
        
        try:
            actual_config = get_diester_config_from_smiles(smiles)
            if actual_config != claimed_config:
                errors.append(
                    f"For option {option}, the provided answer claims the diester is '{claimed_config}', "
                    f"but analysis of its SMILES string reveals it is actually '{actual_config}'."
                )
        except ValueError as e:
            errors.append(f"Error analyzing option {option}: {e}")

    if not errors:
        # This case is unlikely given the known error, but included for completeness.
        # Even if the cis/trans claims were correct, the endo/exo claims would also need verification.
        # However, since the first step of reasoning is flawed, the whole argument fails.
        return "Correct"
    else:
        error_summary = "\n".join(errors)
        reason = (
            "The provided answer is incorrect because its reasoning is based on a flawed premise.\n"
            f"{error_summary}\n"
            "The initial step of the provided answer's logic is to filter for 'cis' isomers, which is chemically correct. "
            "However, it incorrectly classifies option D as 'trans' and eliminates it. "
            "Because this fundamental classification is wrong, the entire logical deduction is invalid, and the final answer 'C' is not reliably supported by the given explanation."
        )
        return reason

# Execute the check and print the result
result = check_answer()
print(result)