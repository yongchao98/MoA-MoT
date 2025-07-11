import sys

# A global flag to stop the search once the first solution is found
solution_found = False

def find_smallest_non_guillotine_tiling():
    """
    This function orchestrates the search for the smallest non-guillotine tiling.
    Based on tiling theory, we target the 5x6 rectangle as the most likely candidate.
    """
    W, H = 5, 6
    # The set of tiles required to tile a 5x6 rectangle (area 30)
    # is two 3x3 squares (area 18) and three 2x2 squares (area 12).
    tiles = [3, 3, 2, 2, 2]
    
    board = [[0] * W for _ in range(H)]
    
    print("Searching for a non-guillotine tiling of a 5x6 rectangle...")
    
    # Start the recursive search for a tiling
    search_tiling_placements(W, H, tiles, board, 1)
    
    if not solution_found:
        print("No non-guillotine tiling was found by this search.")

def search_tiling_placements(W, H, remaining_tiles, board, next_tile_id):
    """
    Performs a backtracking search to find a tiling.
    It works by finding the first empty cell and trying to place each available tile there.
    """
    global solution_found
    if solution_found:
        return

    # If no tiles are left, we have successfully tiled the whole board.
    if not remaining_tiles:
        if not is_guillotine(board, W, H):
            solution_found = True
            print_solution(board, W, H, [3, 3, 2, 2, 2])
        return

    # Find the first empty cell (top-most, then left-most).
    r, c = -1, -1
    for i in range(H):
        for j in range(W):
            if board[i][j] == 0:
                r, c = i, j
                break
        if r != -1:
            break

    # If board is full but we have tiles left, something is wrong. Should not happen if areas match.
    if r == -1:
        return

    # Try placing each unique remaining tile size at the found empty spot.
    tried_sizes = set()
    for i, s in enumerate(remaining_tiles):
        if s in tried_sizes:
            continue
        tried_sizes.add(s)

        if can_place(board, W, H, r, c, s):
            place(board, r, c, s, next_tile_id)
            
            # Recurse with the remaining tiles
            new_remaining_tiles = remaining_tiles[:i] + remaining_tiles[i+1:]
            search_tiling_placements(W, H, new_remaining_tiles, board, next_tile_id + 1)
            
            if solution_found:
                return
            
            # Backtrack: remove the tile by resetting the area to 0
            place(board, r, c, s, 0)

def can_place(board, W, H, r, c, s):
    """Checks if a square of size s can be placed at (r, c)."""
    if r + s > H or c + s > W:
        return False
    for i in range(s):
        for j in range(s):
            if board[r + i][c + j] != 0:
                return False
    return True

def place(board, r, c, s, tile_id):
    """Places or removes a tile on the board."""
    for i in range(s):
        for j in range(s):
            board[r + i][c + j] = tile_id

def is_guillotine(board, W, H):
    """Checks if a tiling is guillotine-cuttable."""
    # Check for any possible vertical cut
    for c in range(1, W):
        is_cut_possible = True
        for r in range(H):
            # If a tile crosses the line between column c-1 and c, the cut is blocked.
            # This is true if the cells on either side belong to the same tile (have same non-zero ID).
            if board[r][c-1] != 0 and board[r][c-1] == board[r][c]:
                is_cut_possible = False
                break
        if is_cut_possible:
            return True

    # Check for any possible horizontal cut
    for r in range(1, H):
        is_cut_possible = True
        for c in range(W):
            if board[r-1][c] != 0 and board[r-1][c] == board[r][c]:
                is_cut_possible = False
                break
        if is_cut_possible:
            return True
            
    return False

def print_solution(board, W, H, tiles_used):
    """Formats and prints the final answer."""
    num_3 = tiles_used.count(3)
    num_2 = tiles_used.count(2)
    area = W * H

    print("\n" + "="*40)
    print("Found a solution!")
    print(f"The smallest rectangle has dimensions {W}x{H}.")
    print(f"The area of this rectangle is {area}.")
    
    print("\nThe area equation, showing each number, is as follows:")
    print(f"{W} * {H} = {area}")
    print(f"{num_3} * 3^2 + {num_2} * 2^2 = {area}")
    print(f"{num_3} * 9 + {num_2} * 4 = {area}")
    print(f"{num_3 * 9} + {num_2 * 4} = {area}")
    print(f"{num_3 * 9 + num_2 * 4} = {area}")

    print("\nA possible non-guillotine tiling is shown below.")
    print("Each number represents a unique square tile:")
    max_id = 0
    for row in board:
        for cell in row:
            if cell > max_id:
                max_id = cell

    width = len(str(max_id))
    for r in range(H):
        print("  ".join(f"{val:>{width}}" for val in board[r]))
    print("="*40)

if __name__ == '__main__':
    find_smallest_non_guillotine_tiling()
    
<<<30>>>