# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def check_the_answer():
    """
    Checks if the proposed answer (4-chlorobenzoic acid) matches all spectral data.
    """
    # 1. Define constraints derived from the spectral data in the question
    constraints = {
        "mw": 156,  # From MS M+ peak m/z = 156
        "has_chlorine": True,  # From MS M+2 peak ratio ~3:1
        "is_carboxylic_acid": True,  # From IR broad O-H and NMR 11.0 ppm peak
        "aromatic_h_count": 4,  # From NMR integration (2H + 2H)
        "aromatic_h_environments": 2,  # From NMR splitting (two doublets) indicating para-substitution
    }

    # 2. Define the proposed answer's structure
    # The final answer from the LLM is D) 4-chlorobenzoic acid
    proposed_answer = {
        "name": "4-chlorobenzoic acid",
        "smiles": "O=C(O)c1ccc(Cl)cc1" 
    }

    # 3. Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(proposed_answer["smiles"])
    if not mol:
        return "Error: Could not parse the SMILES string for the proposed answer."

    # 4. Perform checks against the constraints

    # Check 4a: Mass Spectrometry - Molecular Weight and Chlorine presence
    # Descriptors.ExactMolWt calculates the mass using the most abundant isotopes (e.g., 12C, 1H, 35Cl)
    exact_mw = Descriptors.ExactMolWt(mol)
    if not (constraints["mw"] - 1 < exact_mw < constraints["mw"] + 1):
        return (f"Incorrect: The molecular weight of {proposed_answer['name']} is {exact_mw:.2f}, "
                f"which does not match the M+ peak at m/z={constraints['mw']}.")

    has_cl = any(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    if constraints["has_chlorine"] and not has_cl:
        return (f"Incorrect: The MS data indicates the presence of chlorine, but the structure for "
                f"{proposed_answer['name']} does not contain chlorine.")

    # Check 4b: IR/NMR - Presence of a carboxylic acid functional group
    # Use a SMARTS pattern to find the -COOH group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)
    if constraints["is_carboxylic_acid"] and not has_carboxylic_acid:
        return (f"Incorrect: The IR and NMR data strongly indicate a carboxylic acid, but "
                f"{proposed_answer['name']} does not have this functional group.")

    # Check 4c: NMR - Aromatic proton count and symmetry
    # Add explicit hydrogens to the molecule for accurate counting
    mol_h = Chem.AddHs(mol)
    
    aromatic_h_count = 0
    for atom in mol_h.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetNeighbors()[0].GetIsAromatic():
            aromatic_h_count += 1
    
    if aromatic_h_count != constraints["aromatic_h_count"]:
        return (f"Incorrect: The NMR data shows {constraints['aromatic_h_count']} aromatic protons, "
                f"but {proposed_answer['name']} has {aromatic_h_count}.")

    # Use canonical atom ranking to determine the number of unique aromatic proton environments
    ranks = list(Chem.CanonicalRankAtoms(mol_h, breakTies=True))
    aromatic_h_ranks = set()
    for atom in mol_h.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetNeighbors()[0].GetIsAromatic():
            aromatic_h_ranks.add(ranks[atom.GetIdx()])
    
    num_aromatic_environments = len(aromatic_h_ranks)
    
    if num_aromatic_environments != constraints["aromatic_h_environments"]:
        return (f"Incorrect: The NMR pattern (two doublets) indicates {constraints['aromatic_h_environments']} "
                f"unique aromatic proton environments due to para-symmetry. However, {proposed_answer['name']} "
                f"has {num_aromatic_environments}, which would produce a different splitting pattern.")

    # 5. Final Verdict
    # If all checks are passed, the answer is correct.
    return "Correct"

# Run the check
result = check_the_answer()
print(result)