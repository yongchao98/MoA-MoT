def find_unknot_indices():
    """
    This function returns the list of indices for knots that are isotopic to the unknot.
    
    The analysis of the knots is as follows:
    - K_1: Not the unknot (It's the 6_1 knot).
    - K_2: Isotopic to the unknot. It can be simplified using Reidemeister moves.
    - K_3: Not the unknot (It's the trefoil knot, 3_1).
    - K_4: Not the unknot (It's the figure-eight knot, 4_1).
    - K_5: Isotopic to the unknot. It can be simplified using Reidemeister moves.
    - K_6: Not the unknot (It's the 6_2 knot).
    
    Therefore, the indices of the unknots are 2 and 5.
    """
    
    # The indices of the knots isotopic to the unknot.
    unknot_indices = [2, 5]
    
    print(unknot_indices)

find_unknot_indices()