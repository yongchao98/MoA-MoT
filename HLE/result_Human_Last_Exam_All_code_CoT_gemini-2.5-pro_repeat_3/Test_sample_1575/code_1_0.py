def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given rules.
    """
    num_elements = 100
    num_classes = 5
    elements_per_class = num_elements // num_classes

    # The problem reduces to calculating the number of elements that must cross
    # the boundaries between the 5 equivalence classes of positions.
    # An adjacent swap moves an element across one such boundary.
    #
    # We have 4 boundaries: (0,1), (1,2), (2,3), (3,4).
    # The total moves is the sum of swaps across each boundary.
    #
    # Destination class for an element starting in class k is (4-k).
    #
    # Moves across boundary k <-> k+1:
    # We count how many sets of 20 elements must cross this boundary.
    # An element starting in class c_start crosses if c_start <= k and its
    # destination (4-c_start) is > k.

    moves_per_boundary = []
    total_moves = 0
    
    # Loop through each boundary (k, k+1) from k=0 to 3
    for k in range(num_classes - 1):
        crossings = 0
        # Check for each starting class c_start if it needs to cross boundary k
        for c_start in range(num_classes):
            c_end = 4 - c_start
            # An element crosses from left to right if it starts at or left of k
            # and ends up right of k.
            if c_start <= k and c_end > k:
                crossings += elements_per_class
        moves_per_boundary.append(crossings)
        total_moves += crossings

    # Print the final equation as requested.
    equation_parts = [str(m) for m in moves_per_boundary]
    equation = " + ".join(equation_parts)
    print(f"The total number of moves is the sum of moves across each class boundary.")
    print(f"Moves = {equation} = {total_moves}")

solve_reversal_moves()