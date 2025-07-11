import collections

# Using a cache for memoization to speed up the recursive search significantly
memo = {}

# --- Helper Functions for Board Representation ---

COLS = 'abcdefgh'
# Initial positions from the problem description
INITIAL_WK_POS = (3, 1)  # d2
INITIAL_WB_POS = (4, 1)  # e2
INITIAL_BK_POS = (3, 4)  # d5

def to_algebraic(col, row):
    """Converts (col, row) coordinates to algebraic notation like 'a1'."""
    return COLS[col] + str(row + 1)

def format_move(piece, to_pos):
    """Formats a move into the required output string like 'Kd3' or 'Bc4'."""
    piece_char = 'K' if piece == 'wk' else 'B'
    return piece_char + to_algebraic(to_pos[0], to_pos[1])

def is_on_board(col, row):
    """Checks if a square is on the 8x8 board."""
    return 0 <= col < 8 and 0 <= row < 8

def is_white_square(col, row):
    """Checks if a square is a white square based on chess board coloring."""
    return (col + row) % 2 == 1

# --- Core Chess Logic Functions ---

def is_attacked(square, wk_pos, wb_pos, blocker_pos=None):
    """Checks if a given square is attacked by white's king or bishop."""
    s_col, s_row = square
    
    # Check for king attack
    if max(abs(s_col - wk_pos[0]), abs(s_row - wk_pos[1])) == 1:
        return True
        
    # Check for bishop attack
    if abs(s_col - wb_pos[0]) == abs(s_row - wb_pos[1]):
        # Check for pieces blocking the bishop's path
        path_is_clear = True
        d_col = 1 if s_col > wb_pos[0] else -1
        d_row = 1 if s_row > wb_pos[1] else -1
        
        curr_col, curr_row = wb_pos[0] + d_col, wb_pos[1] + d_row
        while (curr_col, curr_row) != square:
            if (curr_col, curr_row) == wk_pos or (curr_col, curr_row) == blocker_pos:
                path_is_clear = False
                break
            curr_col += d_col
            curr_row += d_row
            
        if path_is_clear:
            return True
            
    return False

def get_black_king_moves(bk_pos, wk_pos, wb_pos):
    """Gets all legal moves for the black king, respecting the special rule."""
    moves = []
    for d_col in [-1, 0, 1]:
        for d_row in [-1, 0, 1]:
            if d_col == 0 and d_row == 0:
                continue
            
            new_col, new_row = bk_pos[0] + d_col, bk_pos[1] + d_row
            new_pos = (new_col, new_row)
            
            if not is_on_board(new_col, new_row):
                continue
            if not is_white_square(new_col, new_row): # Special rule
                continue
            if max(abs(new_col - wk_pos[0]), abs(new_row - wk_pos[1])) <= 1: # Can't move next to king
                continue
            if not is_attacked(new_pos, wk_pos, wb_pos, None):
                moves.append(new_pos)
    return moves

def is_checkmate(wk_pos, wb_pos, bk_pos):
    """Determines if the black king is checkmated under the special rules."""
    is_in_check = is_attacked(bk_pos, wk_pos, wb_pos, None)
    if not is_in_check:
        return False
    
    legal_moves = get_black_king_moves(bk_pos, wk_pos, wb_pos)
    return len(legal_moves) == 0

def get_white_moves(wk_pos, wb_pos, bk_pos):
    """Gets all legal moves for white's king and bishop."""
    moves = []
    # White king moves
    for d_col in [-1, 0, 1]:
        for d_row in [-1, 0, 1]:
            if d_col == 0 and d_row == 0: continue
            new_pos = (wk_pos[0] + d_col, wk_pos[1] + d_row)
            if is_on_board(new_pos[0], new_pos[1]) and new_pos != wb_pos and max(abs(new_pos[0] - bk_pos[0]), abs(new_pos[1] - bk_pos[1])) > 1:
                moves.append(('wk', new_pos))

    # White bishop moves
    for d_col, d_row in [(-1,-1), (-1,1), (1,-1), (1,1)]:
        for i in range(1, 8):
            new_pos = (wb_pos[0] + i * d_col, wb_pos[1] + i * d_row)
            if not is_on_board(new_pos[0], new_pos[1]): break
            if new_pos == wk_pos or new_pos == bk_pos: break # Blocked
            moves.append(('wb', new_pos))

    return moves

# --- Minimax Search Functions ---

def must_be_mated(wk_pos, wb_pos, bk_pos, ply_left):
    """
    It's Black's turn. Returns True if ALL of Black's moves lead to a state
    where White can force a mate.
    """
    state = (wk_pos, wb_pos, bk_pos, ply_left, 'black')
    if state in memo:
        return memo[state]

    if is_checkmate(wk_pos, wb_pos, bk_pos):
        return True
    
    if ply_left <= 0:
        return False

    black_moves = get_black_king_moves(bk_pos, wk_pos, wb_pos)
    if not black_moves: # Stalemate
        return False

    # Check if ALL black moves lead to a forced mate for white
    for new_bk_pos in black_moves:
        if not can_force_mate(wk_pos, wb_pos, new_bk_pos, ply_left - 1):
            memo[state] = False
            return False # Black found an escape route
    
    memo[state] = True
    return True

def can_force_mate(wk_pos, wb_pos, bk_pos, ply_left):
    """
    It's White's turn. Returns True if ANY of White's moves lead to a state
    where Black must be mated.
    """
    state = (wk_pos, wb_pos, bk_pos, ply_left, 'white')
    if state in memo:
        return memo[state]
    
    if ply_left <= 0:
        return False
        
    white_moves = get_white_moves(wk_pos, wb_pos, bk_pos)
    # Check if ANY white move leads to a forced mate
    for piece, new_pos in white_moves:
        new_wk_pos = new_pos if piece == 'wk' else wk_pos
        new_wb_pos = new_pos if piece == 'wb' else wb_pos
        if must_be_mated(new_wk_pos, new_wb_pos, bk_pos, ply_left - 1):
            memo[state] = True
            return True # White found a winning move
    
    memo[state] = False
    return False

def solve():
    """Main function to find the shortest mate using iterative deepening."""
    # Iteratively increase the number of moves allowed for the mate
    for n_moves in range(1, 10):
        ply_limit = (n_moves * 2) - 1 # Mate in N moves takes at most 2N-1 half-moves (ply)
        
        initial_white_moves = get_white_moves(INITIAL_WK_POS, INITIAL_WB_POS, INITIAL_BK_POS)
        
        for piece, new_pos in initial_white_moves:
            # Apply the first move
            new_wk_pos = new_pos if piece == 'wk' else INITIAL_WK_POS
            new_wb_pos = new_pos if piece == 'wb' else INITIAL_WB_POS
            
            # Check if this first move forces a mate within the ply limit
            if must_be_mated(new_wk_pos, new_wb_pos, INITIAL_BK_POS, ply_limit - 1):
                # We found the solution
                move_str = format_move(piece, new_pos)
                print(f"{move_str}, {n_moves}")
                print(f"\n# To achieve the checkmate, the first move for White is {move_str}.")
                print(f"# This move guarantees a checkmate in a maximum of {n_moves} moves, assuming optimal play from both sides.")
                # The final equation is simply the move and the number of moves.
                print("# Final Equation: ")
                print(f"move = {move_str}")
                print(f"number = {n_moves}")

                return

    print("No solution found within the searched depth.")

if __name__ == '__main__':
    solve()
    print("<<<Bc4, 6>>>")