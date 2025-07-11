import sys

def solve_reagents():
    """
    This function analyzes the chemical reaction and identifies the unknown reagents A and B.
    """

    # Analysis for Reagent A
    # The transformation from compound 1 to 2 shows the replacement of a ring oxygen atom
    # with a H2N-N< group. This is a characteristic reaction of pyrylium-type salts with hydrazine.
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    # Analysis for Reagent B
    # The transformation from compound 2 to 3 involves two main changes:
    # 1. The N-NH2 group is converted to an NH group.
    # 2. A methoxy group is replaced by an N-propyl group.
    # The source for the N-propyl group must be reagent B.
    reagent_B = {
        "name": "n-Propylamine",
        "formula": "CH3CH2CH2NH2"
    }

    # Print the identified reagents
    print(f"Based on the reaction scheme analysis:")
    print(f"Reagent A is: {reagent_A['name']} ({reagent_A['formula']})")
    print(f"Reagent B is: {reagent_B['name']} ({reagent_B['formula']})")

solve_reagents()