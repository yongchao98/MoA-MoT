def find_unknot_indices():
    """
    Analyzes six knot diagrams and identifies which ones are isotopic to the unknot.

    The analysis for each knot is as follows:
    - K1: A non-trivial knot (8_18 knot). Not an unknot.
    - K2: A known complex diagram of the unknot (a "hard unknot"). It is an unknot.
    - K3: A simple diagram that simplifies to the unknot via one Reidemeister II move. It is an unknot.
    - K4: The figure-eight knot (4_1). Not an unknot.
    - K5: Simplifies to the unknot via two sequential Reidemeister II moves. It is an unknot.
    - K6: A non-trivial knot (6_2 knot). Not an unknot.
    """
    
    # A list representing whether each knot K_i is an unknot.
    # Indices correspond to K1, K2, ..., K6.
    is_unknot = [False, True, True, False, True, False]
    
    unknot_indices = []
    for i in range(len(is_unknot)):
        if is_unknot[i]:
            # The knots are 1-indexed, so we add 1 to the loop index.
            unknot_indices.append(i + 1)
            
    print(unknot_indices)

find_unknot_indices()