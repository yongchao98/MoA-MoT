from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import re

def find_heavier_product():
    """
    Identifies the product with the higher molar mass from the hydrolysis of a plausible ketal.
    The ketal is assumed to be formed from acetophenone and 1,1-bis(hydroxymethyl)cyclobutane,
    as inferred from the malformed SMILES string and reaction conditions provided in the prompt.
    """
    # Define the plausible products of the hydrolysis reaction
    product1_smi = 'CC(=O)c1ccccc1'
    product1_name = 'Acetophenone'
    product2_smi = 'C1CC1(CO)CO'
    product2_name = '1,1-bis(hydroxymethyl)cyclobutane'

    products = [(product1_name, product1_smi), (product2_name, product2_smi)]
    masses = []

    # Use atomic weights of the most common isotopes for exact mass calculation
    atomic_weights = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
    }

    print("The reaction is assumed to be a ketal hydrolysis.")
    print("The predicted products are Acetophenone and 1,1-bis(hydroxymethyl)cyclobutane.")
    print("\nCalculating the molar mass for each product...")
    print("-" * 60)

    for name, smi in products:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Error: Could not parse SMILES for {name}: {smi}")
            continue

        # Get the molecular formula from RDKit
        formula = rdMolDescriptors.CalcMolFormula(mol)
        
        # Calculate the exact molar weight using RDKit
        exact_mass = Descriptors.ExactMolWt(mol)
        masses.append(exact_mass)
        
        print(f"Product: {name}")
        print(f"SMILES: {smi}")
        print(f"Formula: {formula}")

        # Parse the formula to get atom counts for the calculation explanation
        atom_counts = {el: int(co) if co else 1 for el, co in re.findall(r'([A-Z][a-z]?)(\d*)', formula)}
        
        calc_str_parts = []
        manual_mass = 0.0
        
        # Build the string showing the calculation for each atom
        for atom, count in atom_counts.items():
            weight = atomic_weights.get(atom, 0)
            mass_contribution = count * weight
            manual_mass += mass_contribution
            calc_str_parts.append(f"({count} * {weight:.6f})") # Show numbers for the equation
        
        calc_equation = f"Equation: {' + '.join(calc_str_parts)} = {manual_mass:.6f}"
        
        print(calc_equation)
        print(f"Exact Molar Mass: {exact_mass:.6f}")
        print("-" * 60)

    # Compare the calculated masses and determine the product with the higher molar mass
    if not masses or len(masses) != 2:
        print("Could not complete the comparison due to an error.")
        return

    if masses[0] > masses[1]:
        winner_smi = products[0][1]
        winner_name = products[0][0]
    else:
        winner_smi = products[1][1]
        winner_name = products[1][0]

    print(f"\nComparison: {masses[0]:.4f} ({products[0][0]}) vs. {masses[1]:.4f} ({products[1][0]})")
    print(f"The product with the higher molar mass is {winner_name}.")
    print("\nSMILES string of the heavier product:")
    print(winner_smi)


if __name__ == '__main__':
    find_heavier_product()
