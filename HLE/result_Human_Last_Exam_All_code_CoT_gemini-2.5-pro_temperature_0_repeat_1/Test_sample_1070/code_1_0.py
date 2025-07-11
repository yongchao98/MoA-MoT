import sys

# Increase recursion limit for deep searches, just in case.
sys.setrecursionlimit(2000)

# --- Helper Functions ---

def to_coords(alg):
    """Converts algebraic notation (e.g., 'd5') to coordinates (3, 4)."""
    col = ord(alg[0]) - ord('a')
    row = int(alg[1]) - 1
    return (col, row)

def to_alg(coords):
    """Converts coordinates (3, 4) to algebraic notation 'd5'."""
    col, row = coords
    return chr(ord('a') + col) + str(row + 1)

def is_white_square(coords):
    """Checks if a square is white. A square (col, row) is white if (col + row) is odd."""
    col, row = coords
    return (col + row) % 2 != 0

def is_attacked(sq_coords, wk_coords, wb_coords):
    """Checks if a square is attacked by the White King or Bishop."""
    # Attacked by King
    if max(abs(sq_coords[0] - wk_coords[0]), abs(sq_coords[1] - wk_coords[1])) == 1:
        return True
    # Attacked by Bishop
    if abs(sq_coords[0] - wb_coords[0]) == abs(sq_coords[1] - wb_coords[1]):
        return True
    return False

# --- Move Generation ---

def generate_white_moves(wk_pos, wb_pos, bk_pos):
    """Generates all legal moves for White."""
    moves = {}
    # King moves
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            nx, ny = wk_pos[0] + dx, wk_pos[1] + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                new_pos = (nx, ny)
                # King cannot move to an occupied square or adjacent to the other king
                if new_pos != wb_pos and max(abs(new_pos[0] - bk_pos[0]), abs(new_pos[1] - bk_pos[1])) > 1:
                    moves[f"K{to_alg(new_pos)}"] = (new_pos, wb_pos)
    # Bishop moves
    for dx, dy in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        for i in range(1, 8):
            nx, ny = wb_pos[0] + i * dx, wb_pos[1] + i * dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                new_pos = (nx, ny)
                # Path is blocked by own king
                if new_pos == wk_pos:
                    break
                # Path is blocked by enemy king
                path_obstructed = False
                for j in range(1, i):
                    path_x, path_y = wb_pos[0] + j * dx, wb_pos[1] + j * dy
                    if (path_x, path_y) == bk_pos:
                        path_obstructed = True
                        break
                if path_obstructed:
                    break
                moves[f"B{to_alg(new_pos)}"] = (wk_pos, new_pos)
            else:
                break
    return moves

def generate_black_moves(wk_pos, wb_pos, bk_pos):
    """Generates all legal moves for the Black King."""
    moves = {}
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            nx, ny = bk_pos[0] + dx, bk_pos[1] + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                new_pos = (nx, ny)
                # Rule: Must move to a white square
                if not is_white_square(new_pos):
                    continue
                # Cannot move to an occupied square
                if new_pos == wk_pos or new_pos == wb_pos:
                    continue
                # Cannot move into an attacked square (check)
                if is_attacked(new_pos, wk_pos, wb_pos):
                    continue
                moves[f"K{to_alg(new_pos)}"] = new_pos
    return moves

# --- Minimax Solver with Memoization ---

memo = {}

def solve_plies(wk_pos, wb_pos, bk_pos, turn):
    """
    Recursively finds the number of plies to a forced mate.
    - White (turn=0) wants to minimize this number.
    - Black (turn=1) wants to maximize this number.
    Returns float('inf') for a draw or an unproven line.
    """
    state = (wk_pos, wb_pos, bk_pos, turn)
    if state in memo:
        return memo[state]

    # White's turn (minimizing player)
    if turn == 0:
        min_plies = float('inf')
        white_moves = generate_white_moves(wk_pos, wb_pos, bk_pos)
        if not white_moves: # Should not happen unless white is stalemated
            memo[state] = float('inf')
            return float('inf')

        for move, (new_wk, new_wb) in white_moves.items():
            plies = solve_plies(new_wk, new_wb, bk_pos, 1)
            min_plies = min(min_plies, plies)
        
        result = float('inf') if min_plies == float('inf') else 1 + min_plies
        memo[state] = result
        return result

    # Black's turn (maximizing player)
    else: # turn == 1
        max_plies = -1
        black_moves = generate_black_moves(wk_pos, wb_pos, bk_pos)
        
        if not black_moves:
            if is_attacked(bk_pos, wk_pos, wb_pos):
                # Checkmate: 0 plies from this point
                memo[state] = 0
                return 0
            else:
                # Stalemate: Draw, infinite plies
                memo[state] = float('inf')
                return float('inf')

        for move, new_bk in black_moves.items():
            plies = solve_plies(wk_pos, wb_pos, new_bk, 0)
            if plies != float('inf'):
                max_plies = max(max_plies, plies)

        if max_plies == -1: # Black can always escape to a line with no forced mate
            result = float('inf')
        else:
            result = 1 + max_plies
            
        memo[state] = result
        return result

# --- Main Execution ---

def find_best_move_and_mate():
    """
    Finds the best first move for White and the number of moves in the forced mate.
    """
    initial_wk = to_coords('d2')
    initial_wb = to_coords('e2')
    initial_bk = to_coords('d5')

    best_move = None
    min_total_plies = float('inf')

    # Iterate through all possible first moves for White
    white_moves = generate_white_moves(initial_wk, initial_wb, initial_bk)
    for move, (new_wk, new_wb) in white_moves.items():
        # After White's move, it's Black's turn. Find how long Black can delay the mate.
        plies_from_black = solve_plies(new_wk, new_wb, initial_bk, 1)
        
        if plies_from_black != float('inf'):
            total_plies = 1 + plies_from_black
            if total_plies < min_total_plies:
                min_total_plies = total_plies
                best_move = move

    if best_move is None:
        print("No forced mate found.")
    else:
        # Convert total plies to number of White's moves
        num_white_moves = (min_total_plies + 1) // 2
        print(f"{best_move}, {num_white_moves}")

if __name__ == '__main__':
    find_best_move_and_mate()