def find_unknot_indices():
    """
    This function returns the indices of the knots that are isotopic to the unknot.

    The analysis is done visually based on knot theory principles:
    - K1: Stevedore knot (6_1), not the unknot.
    - K2: A complex diagram that simplifies to the unknot.
    - K3: A diagram that simplifies to the unknot after pulling out the bottom strand.
    - K4: Figure-eight knot (4_1), not the unknot.
    - K5: A diagram that simplifies to the unknot after two Reidemeister II moves.
    - K6: Cinquefoil knot (5_1), not the unknot.
    """
    unknot_indices = [2, 3, 5]
    print(unknot_indices)

find_unknot_indices()