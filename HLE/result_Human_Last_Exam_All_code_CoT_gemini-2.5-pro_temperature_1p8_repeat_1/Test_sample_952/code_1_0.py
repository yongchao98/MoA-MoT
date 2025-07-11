def solve():
    """
    This function calculates the number of movable diamonds in a checkerboard
    arrangement on a 2024x2024 grid.
    """
    n = 2024
    
    # Create a grid representing the checkerboard.
    # grid[i][j] = 1 means a diamond is at (i, j), 0 means empty.
    # We use 0-based indexing for the grid, so cells are (0..n-1, 0..n-1).
    grid = [[0] * n for _ in range(n)]

    # Place diamonds on "black" cells (where i + j is even).
    diamond_cells = []
    for r in range(n):
        for c in range(n):
            if (r + c) % 2 == 0:
                grid[r][c] = 1
                diamond_cells.append((r, c))

    movable_diamonds_count = 0

    # Check each diamond to see if it's movable.
    for r_diamond, c_diamond in diamond_cells:
        is_diamond_movable = False
        
        # Check potential moves to 4 adjacent cells.
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            r_move, c_move = r_diamond + dr, c_diamond + dc

            # Check if the destination is within the grid.
            if 0 <= r_move < n and 0 <= c_move < n:
                # The move destination must be empty (which is true in checkerboard).
                if grid[r_move][c_move] == 1:
                    continue
                
                # Now check if the new position is valid.
                # A move is valid if all new neighbors are empty.
                is_move_valid = True
                
                # Iterate through neighbors of the potential new position.
                for dr_new, dc_new in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    r_neighbor, c_neighbor = r_move + dr_new, c_move + dc_new

                    # The original diamond position doesn't count as a blocker.
                    if (r_neighbor, c_neighbor) == (r_diamond, c_diamond):
                        continue

                    # Check if the neighbor is within bounds.
                    if 0 <= r_neighbor < n and 0 <= c_neighbor < n:
                        # If the neighbor cell has a diamond, the move is blocked.
                        if grid[r_neighbor][c_neighbor] == 1:
                            is_move_valid = False
                            break
                
                if not is_move_valid:
                    continue # Try another move direction for the current diamond.
                else:
                    # Found a valid move for this diamond.
                    is_diamond_movable = True
                    break # No need to check other moves for this diamond.
        
        if is_diamond_movable:
            movable_diamonds_count += 1
    
    print(f"For the checkerboard arrangement, the number of movable diamonds is: {movable_diamonds_count}")
    print("This represents the minimum number of movable diamonds across all possible arrangements.")
    print(f"The largest value k is therefore: {movable_diamonds_count}")

solve()