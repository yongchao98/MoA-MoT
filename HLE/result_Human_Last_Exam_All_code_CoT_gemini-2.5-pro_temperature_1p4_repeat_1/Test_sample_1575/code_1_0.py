def solve_reversal_moves():
    """
    Calculates the minimum number of adjacent swaps to reverse a sequence of 100 elements
    with free non-adjacent swaps at a distance of 5.
    """
    N = 100
    group_count = 5
    group_size = N // group_count

    print(f"The sequence has {N} elements.")
    print(f"Free swaps of distance 5 partition the positions into {group_count} groups of {group_size} elements each.")
    print("Let's analyze the required movement of elements between these groups.")
    print("\nAn element at position i must move to 101 - i.")
    print("This defines a mapping between the initial group of an element and its target group.")
    print("The mapping is G_k -> G_(4-k) for groups k = 0, 1, 2, 3, 4.")
    print("This results in the following exchanges:")
    print("- Elements from Group 0 must exchange with elements from Group 4.")
    print("- Elements from Group 1 must exchange with elements from Group 3.")
    print("- Elements from Group 2 stay within Group 2 (cost = 0).\n")

    # Cost for G0 <-> G4 exchange
    # This crosses one boundary (between G0 and G4)
    # Flow is 20 in each direction. Moves = max(20, 20).
    moves_g0_g4 = group_size
    print("Cost for G0 <-> G4 exchange:")
    print(f"This requires moving {group_size} elements across one boundary.")
    print(f"Moves = {moves_g0_g4}\n")

    # Cost for G1 <-> G3 exchange
    # This crosses two boundaries (G1-G2 and G2-G3)
    # At boundary G1-G2, flow is 20 each way.
    moves_g1_g2_boundary = group_size
    # At boundary G2-G3, flow is 20 each way.
    moves_g2_g3_boundary = group_size
    total_moves_g1_g3 = moves_g1_g2_boundary + moves_g2_g3_boundary
    print("Cost for G1 <-> G3 exchange:")
    print(f"This requires crossing two boundaries (G1-G2 and G2-G3).")
    print(f"Moves at G1-G2 boundary = {moves_g1_g2_boundary}")
    print(f"Moves at G2-G3 boundary = {moves_g2_g3_boundary}")
    print(f"Total for this exchange = {total_moves_g1_g3}\n")

    # Total moves
    total_moves = moves_g0_g4 + total_moves_g1_g3

    print("The total minimum number of moves is the sum of these costs.")
    print(f"Total moves = {moves_g0_g4} (for G0<->G4) + {moves_g1_g2_boundary} (at G1-G2) + {moves_g2_g3_boundary} (at G2-G3)")
    print(f"Final calculation: {moves_g0_g4} + {moves_g1_g2_boundary} + {moves_g2_g3_boundary} = {total_moves}")

solve_reversal_moves()