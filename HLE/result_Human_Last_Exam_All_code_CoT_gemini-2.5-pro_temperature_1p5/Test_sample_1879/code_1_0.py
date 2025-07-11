def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint that no more than 3 consecutive steps in the same
    direction are allowed.
    """
    R_MAX = 4
    U_MAX = 8
    CONSECUTIVE_LIMIT = 3

    # H[r][u]: ways to get to (r, u) ending with a Right move
    # V[r][u]: ways to get to (r, u) ending with an Up move
    H = [[0 for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]
    V = [[0 for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]

    # Base cases for axis u=0 (only Right moves possible: R, RR, RRR)
    for r in range(1, R_MAX + 1):
        if r <= CONSECUTIVE_LIMIT:
            H[r][0] = 1

    # Base cases for axis r=0 (only Up moves possible: U, UU, UUU)
    for u in range(1, U_MAX + 1):
        if u <= CONSECUTIVE_LIMIT:
            V[0][u] = 1

    # Fill the DP tables using the recurrence relations
    for r in range(1, R_MAX + 1):
        for u in range(1, U_MAX + 1):
            # To calculate H[r][u], sum the ways to arrive from the left
            # from a position that was reached by an Up move.
            # This can happen by appending 1 to 3 'R's to a path ending in 'U'.
            # e.g., ...UR, ...URR, ...URRR
            for k in range(1, CONSECUTIVE_LIMIT + 1):
                if r - k >= 0:
                    H[r][u] += V[r - k][u]

            # To calculate V[r][u], sum the ways to arrive from below
            # from a position that was reached by a Right move.
            # This can happen by appending 1 to 3 'U's to a path ending in 'R'.
            # e.g., ...RU, ...RUU, ...RUUU
            for k in range(1, CONSECUTIVE_LIMIT + 1):
                if u - k >= 0:
                    V[r][u] += H[r][u - k]

    final_H = H[R_MAX][U_MAX]
    final_V = V[R_MAX][U_MAX]
    total_ways = final_H + final_V

    print(f"To reach ({R_MAX},{U_MAX}):")
    print(f"Number of paths ending with a Right move: {final_H}")
    print(f"Number of paths ending with an Up move: {final_V}")
    print(f"Total unique ways = {final_H} + {final_V} = {total_ways}")

solve_grid_path()