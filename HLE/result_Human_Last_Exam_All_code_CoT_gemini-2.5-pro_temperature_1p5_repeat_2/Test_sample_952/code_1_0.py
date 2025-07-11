def solve():
    """
    This function calculates the minimum number of movable diamonds for a specific
    "checkerboard" arrangement on an N x N grid.
    If we can find an arrangement with 0 movable diamonds, the value k must be 0.
    """
    N = 2024

    # Create the checkerboard arrangement. A cell (r, c) has a diamond if r+c is even.
    diamonds = set()
    for r in range(N):
        for c in range(N):
            if (r + c) % 2 == 0:
                diamonds.add((r, c))

    movable_diamonds_coords = set()

    # Iterate through each diamond to check if it's movable
    for r_d, c_d in diamonds:
        is_diamond_movable = False
        # A diamond can have up to 4 neighbors
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            r_n, c_n = r_d + dr, c_d + dc

            # Check if neighbor is within the grid
            if 0 <= r_n < N and 0 <= c_n < N:
                # The neighbor cell must be empty (which is true in checkerboard)
                if (r_n, c_n) in diamonds:
                    continue

                # Check if moving to this neighbor cell is a valid move.
                # All other neighbors of (r_n, c_n) must be empty.
                is_move_valid = True
                # A neighbor cell can have up to 4 neighbors
                for dr_nn, dc_nn in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    r_nn, c_nn = r_n + dr_nn, c_n + dc_nn

                    # This is the original position of our diamond, which will be empty after the move
                    if (r_nn, c_nn) == (r_d, c_d):
                        continue
                    
                    # Check if neighbor's neighbor is within the grid
                    if 0 <= r_nn < N and 0 <= c_nn < N:
                        # If any other neighbor has a diamond, the move is invalid
                        if (r_nn, c_nn) in diamonds:
                            is_move_valid = False
                            break
                
                if not is_move_valid:
                    continue
                else:
                    # Found a valid move for this diamond
                    is_diamond_movable = True
                    break # No need to check other moves for this diamond
        
        if is_diamond_movable:
            movable_diamonds_coords.add((r_d, c_d))

    k = len(movable_diamonds_coords)
    
    print("For the checkerboard arrangement, the number of movable diamonds is:")
    print(f"M = {k}")
    print("\nSince there exists a valid arrangement with 0 movable diamonds,")
    print("the largest value k such that *every* arrangement has at least k movable diamonds must be 0.")
    print("\nFinal Equation:")
    print(f"k = {k}")

solve()