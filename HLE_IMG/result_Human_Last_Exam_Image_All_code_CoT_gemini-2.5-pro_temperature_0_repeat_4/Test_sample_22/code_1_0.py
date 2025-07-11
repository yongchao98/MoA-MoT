def find_unknot_indices():
    """
    Analyzes the six provided knot diagrams and identifies which are isotopic to the unknot.

    The analysis is as follows:
    - K_1: Stevedore knot (6_1). Not the unknot.
    - K_2: A complex diagram that simplifies to the unknot via Reidemeister moves. It is the unknot.
    - K_3: A diagram that simplifies to the unknot via a single Reidemeister II move. It is the unknot.
    - K_4: Figure-eight knot (4_1). Not the unknot.
    - K_5: A diagram that simplifies to the unknot via Reidemeister II and I moves. It is the unknot.
    - K_6: Trefoil knot (3_1). Not the unknot.

    The indices of the unknots are 2, 3, and 5.
    """
    unknot_indices = [2, 3, 5]
    print(unknot_indices)

find_unknot_indices()