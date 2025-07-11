# First, please make sure you have RDKit and cirpy installed:
# pip install rdkit cirpy

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import cirpy
import re

def solve_reaction():
    """
    Solves the anionic oxy-Cope rearrangement problem for the given substrate.
    """
    reactant_name = '(1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol'

    print(f"Starting Material: {reactant_name}\n")

    # Step 1: Use a web service to convert the chemical name to a SMILES string
    # Using a known SMILES string for reliability as web services can be unstable.
    # smi_reactant = cirpy.resolve(reactant_name, 'smiles')
    smi_reactant = 'CC(C)(C)[Si](C)(C)O[C@H]1C/C=C(\\C1)C2(O)[C@@H]3C=C[C@H]([C@@H]2C3(OC)OC)C'
    
    if not smi_reactant:
        print("Could not resolve the chemical name to a structure.")
        return

    mol_reactant = Chem.MolFromSmiles(smi_reactant)
    if not mol_reactant:
        print("Could not create molecule from SMILES string.")
        return
        
    # The reaction is an isomerization, so the formula remains the same.
    formula = CalcMolFormula(mol_reactant)
    
    print("This reaction is an Anionic Oxy-Cope Rearrangement, which is an isomerization.")
    print("Therefore, the molecular formula of the product is the same as the reactant.\n")
    
    # Define the reaction using a reaction SMARTS (SMIRKS) pattern
    # Reactant pattern: [C:1]=[C:2]-[C:3]([O:7][H])-[C:4]-[C:5]=[C:6]
    # Product pattern after tautomerization: a new C1-C6 bond, C3 becomes C=O, and pi-bonds shift.
    # The SMARTS [C:1]1...[C:6]-1 indicates the new bond between atoms mapped to 1 and 6.
    smirks = '[C:1]=[C:2]-[C:3](-[O:7][H])-[C:4]-[C:5]=[C:6]>>[C:1]1-[C:2]-[C:3](=[O:7])-[C:4]=[C:5]-[C:6]1'
    rxn = AllChem.ReactionFromSmarts(smirks)

    # Run the reaction
    products = rxn.RunReactants((mol_reactant,))

    if not products:
        print("Reaction could not be applied to the starting material.")
        return
        
    # We expect one major product from this intramolecular reaction
    mol_product = products[0][0]
    
    try:
        Chem.SanitizeMol(mol_product)
    except Exception as e:
        print(f"Could not sanitize the product molecule: {e}")
        return

    # Print the reaction equation with molecular formulas
    print("Reaction Equation (by formula):")
    print(f"{formula} --> {CalcMolFormula(mol_product)}\n")

    # Print the numbers for each element as requested
    print("Component atoms and their counts:")
    # Use regex to find element and count pairs
    atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    for element, count in atom_counts:
        # If count is not specified, it's 1
        num = int(count) if count else 1
        print(f"{element}: {num}")
    
    # For chemical context, let's also print the SMILES of the product.
    # We remove stereochemistry because predicting it accurately requires complex transition state analysis.
    Chem.RemoveStereochemistry(mol_product)
    smi_product = Chem.MolToSmiles(mol_product)
    print(f"\nProduct Structure (SMILES, stereochemistry omitted):\n{smi_product}")


solve_reaction()
# The final answer is the chemical formula of the product.
# Since this is an isomerization reaction, the formula doesn't change.
# Let's get the formula from RDKit.
reactant_smi = 'CC(C)(C)[Si](C)(C)O[C@H]1C/C=C(\\C1)C2(O)[C@@H]3C=C[C@H]([C@@H]2C3(OC)OC)C'
reactant_mol = Chem.MolFromSmiles(reactant_smi)
final_formula = CalcMolFormula(reactant_mol)
print(f'<<<{final_formula}>>>')
