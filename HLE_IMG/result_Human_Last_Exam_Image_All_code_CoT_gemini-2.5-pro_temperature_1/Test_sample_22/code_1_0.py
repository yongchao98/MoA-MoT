def find_unknots():
    """
    This function returns the indices of the knots that are isotopic to the unknot.
    The analysis is based on visual inspection and knowledge of knot theory.

    - K1: 6_1 knot (Stevedore knot). Not an unknot.
    - K2: A complex diagram of the unknot, can be "unwoven". It is an unknot.
    - K3: A 2-crossing diagram. All 1- or 2-crossing diagrams are unknots. It is an unknot.
    - K4: 6_2 knot. Not an unknot.
    - K5: 5_2 knot. Not an unknot.
    - K6: Another diagram of the 6_1 knot. Not an unknot.
    """
    unknot_indices = [2, 3]
    print(unknot_indices)

find_unknots()