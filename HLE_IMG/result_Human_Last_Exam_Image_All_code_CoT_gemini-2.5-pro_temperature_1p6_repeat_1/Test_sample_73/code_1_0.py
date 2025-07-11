def solve_stereochemistry():
    """
    Determines and prints the stereochemical assignments for the reaction.
    The order is: Reactant 1, Reactant 2, Product-left-center, Product-right-center.
    """
    # Stereochemical assignment for the first reactant (acid chloride)
    reactant1_config = "(S)"

    # Stereochemical assignment for the second reactant (alcohol)
    reactant2_config = "(R)"

    # The reaction is an esterification that retains the stereochemistry of both centers.
    # The product is drawn with the alcohol-derived part on the left and the acid-derived part on the right.
    product_config_left = reactant2_config  # from alcohol
    product_config_right = reactant1_config # from acid chloride

    # The four assignments moving from left to right in the scheme are:
    # 1. Reactant 1
    # 2. Reactant 2
    # 3. Product's left stereocenter (derived from alcohol)
    # 4. Product's right stereocenter (derived from acid chloride)
    assignments = [
        reactant1_config,
        reactant2_config,
        product_config_left,
        product_config_right,
    ]

    print("The stereochemical assignments for the four stereocenters from left to right are:")
    print(", ".join(assignments))

solve_stereochemistry()