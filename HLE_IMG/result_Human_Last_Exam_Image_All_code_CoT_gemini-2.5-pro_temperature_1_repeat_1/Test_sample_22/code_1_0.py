def find_unknot_indices():
    """
    This function returns the list of indices of knots that are isotopic to the unknot,
    based on a visual analysis of the provided image.

    The analysis is as follows:
    - K1: Not the unknot (it's the 6_1 knot).
    - K2: Is the unknot (a complex diagram that can be simplified).
    - K3: Is the unknot (a 3-crossing diagram of the unknot).
    - K4: Not the unknot (it's the figure-eight knot, 4_1).
    - K5: Is the unknot (simplifies to K3 via a Reidemeister II move).
    - K6: Not the unknot (it's the 5_1 knot).
    """
    unknot_indices = [2, 3, 5]
    print(unknot_indices)

find_unknot_indices()