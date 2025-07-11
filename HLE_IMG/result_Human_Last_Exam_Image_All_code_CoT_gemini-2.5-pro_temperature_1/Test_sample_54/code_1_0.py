def identify_reagents():
    """
    This function identifies the reagents A and B in the provided chemical reaction scheme.
    """
    # Step 1: Identification of Reagent A
    # The transformation from compound 1 to 2 involves the replacement of a ring oxygen atom
    # with an N-NH2 group. This is a characteristic reaction of pyrylium-like salts with hydrazine.
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    # Step 2: Identification of Reagent B
    # The transformation from compound 2 to 3 involves two changes:
    # 1. Replacement of another ring oxygen with an N-propyl group.
    # 2. Conversion of the N-NH2 group to an NH group.
    # These reactions are consistent with the use of n-propylamine, which provides the N-propyl group
    # and facilitates the formation of the quinacridinium core structure.
    reagent_B = {
        "name": "n-Propylamine",
        "formula": "CH3CH2CH2NH2"
    }

    # Output the findings
    print("Based on the reaction scheme, the reagents are identified as follows:")
    print(f"Reagent A is: {reagent_A['name']} (Formula: {reagent_A['formula']})")
    print(f"Reagent B is: {reagent_B['name']} (Formula: {reagent_B['formula']})")

    print("\nThe reaction equations with the identified reagents are:")
    print(f"Equation 1: Compound 1 + {reagent_A['name']} -> Compound 2")
    print(f"Equation 2: Compound 2 + {reagent_B['name']} -> Compound 3")

if __name__ == "__main__":
    identify_reagents()