import requests
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_smiles_from_name(name: str) -> str | None:
    """Converts an IUPAC name to a SMILES string using the NCI/CADD resolver."""
    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        smiles = response.text.strip()
        if "html" in smiles:  # Error page returns html
            return None
        return smiles
    except requests.exceptions.RequestException:
        return None

def get_mol_formula(smiles: str) -> str | None:
    """Calculates the molecular formula from a SMILES string."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)

def find_diene_system(mol, connecting_bonds):
    """Finds a diene system with a specific number of connecting single bonds."""
    num_atoms_in_path = connecting_bonds + 4
    paths = Chem.FindAllPathsOfLengthN(mol, num_atoms_in_path - 1)
    
    for path in paths:
        bond1 = mol.GetBondBetweenAtoms(path[0], path[1])
        bondN = mol.GetBondBetweenAtoms(path[-2], path[-1])
        
        if not (bond1 and bond1.GetBondType() == Chem.BondType.DOUBLE and
                bondN and bondN.GetBondType() == Chem.BondType.DOUBLE):
            continue
            
        is_single = True
        for i in range(1, num_atoms_in_path - 2):
            bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
            if not (bond and bond.GetBondType() == Chem.BondType.SINGLE):
                is_single = False
                break
        
        if is_single:
            return path
            
    return None

def perform_3_3_shift(mol, system_indices):
    """Performs a [3,3]-sigmatropic shift on a 1,7-diene system."""
    # System is C1=C2-C3-C4-C5-C6=C7
    # The rearranging fragments are C2,C3,C4 and C5,C6,C7
    # Break C4-C5, form C2-C7, shift pi bonds to C3=C4 and C5=C6
    c1, c2, c3, c4, c5, c6, c7 = system_indices
    
    rw_mol = Chem.RWMol(mol)
    
    rw_mol.RemoveBond(c4, c5)
    rw_mol.AddBond(c2, c7, Chem.BondType.SINGLE)
    rw_mol.GetBondBetweenAtoms(c2, c3).SetBondType(Chem.BondType.SINGLE)
    rw_mol.GetBondBetweenAtoms(c6, c7).SetBondType(Chem.BondType.SINGLE)
    rw_mol.AddBond(c3, c4, Chem.BondType.DOUBLE)
    rw_mol.AddBond(c5, c6, Chem.BondType.DOUBLE)
    
    product_mol = rw_mol.GetMol()
    
    try:
        Chem.SanitizeMol(product_mol)
        return product_mol, None
    except Exception as e:
        return None, f"Failed to sanitize the product molecule: {e}"

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the chemical reaction.
    """
    start_material_name = "5-butylnona-2,6-diene"
    llm_answer_key = "D"
    llm_answer_name = "4-ethyl-3-methyldeca-1,5-diene"

    # Step 1: Resolve names and check for isomerism
    start_smiles = get_smiles_from_name(start_material_name)
    if not start_smiles:
        return f"Error: Could not resolve the starting material '{start_material_name}' into a chemical structure."

    llm_answer_smiles = get_smiles_from_name(llm_answer_name)
    if not llm_answer_smiles:
        return f"Error: Could not resolve the answer '{llm_answer_name}' into a chemical structure."

    start_formula = get_mol_formula(start_smiles)
    llm_answer_formula = get_mol_formula(llm_answer_smiles)

    if start_formula != llm_answer_formula:
        return (f"Incorrect: The proposed answer (Option {llm_answer_key}) is not an isomer of the starting material.\n"
                f"Starting material formula: {start_formula}\n"
                f"Option {llm_answer_key} formula: {llm_answer_formula}")

    # Step 2: Simulate the reaction
    start_mol = Chem.MolFromSmiles(start_smiles)
    
    # The name "nona-2,6-diene" implies a 1,7-diene system (3 single bonds between double bonds)
    system = find_diene_system(start_mol, connecting_bonds=3)
    
    if not system:
        return "Error: Could not programmatically identify the 1,7-diene system in the starting material."
        
    predicted_mol, error = perform_3_3_shift(start_mol, system)
    
    if error:
        return f"Error during reaction simulation: {error}"

    # Step 3: Compare the predicted product with the LLM's answer
    # Use canonical SMILES without stereochemistry for a robust comparison of connectivity
    predicted_canon_smiles = Chem.MolToSmiles(predicted_mol, isomericSmiles=False)
    llm_answer_canon_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(llm_answer_smiles), isomericSmiles=False)

    if predicted_canon_smiles == llm_answer_canon_smiles:
        return "Correct"
    else:
        return (f"Incorrect: The product derived from the [3,3]-sigmatropic shift mechanism does not match the proposed answer (Option {llm_answer_key}).\n"
                f"  - Predicted product SMILES: {predicted_canon_smiles}\n"
                f"  - Proposed answer SMILES:   {llm_answer_canon_smiles}")

# Run the check
result = check_correctness()
print(result)