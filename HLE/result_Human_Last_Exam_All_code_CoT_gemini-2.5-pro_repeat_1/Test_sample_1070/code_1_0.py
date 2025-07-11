import collections

# Helper functions for chess coordinates
def alg_to_coord(alg):
    """Converts algebraic notation (e.g., 'd5') to (col, row) tuple."""
    col = ord(alg[0]) - ord('a')
    row = int(alg[1]) - 1
    return (col, row)

def coord_to_alg(coord):
    """Converts (col, row) tuple to algebraic notation."""
    col, row = coord
    return chr(ord('a') + col) + str(row + 1)

def piece_to_alg(piece_char, coord):
    """Formats a move string like 'Bd3'."""
    return piece_char + coord_to_alg(coord)

# Board and piece logic
def is_white_square(coord):
    """Checks if a square is white. a1=(0,0) is dark."""
    return (coord[0] + coord[1]) % 2 != 0

def get_king_moves(pos):
    """Generates all potential king moves from a position."""
    col, row = pos
    for dc in [-1, 0, 1]:
        for dr in [-1, 0, 1]:
            if dc == 0 and dr == 0:
                continue
            new_col, new_row = col + dc, row + dr
            if 0 <= new_col < 8 and 0 <= new_row < 8:
                yield (new_col, new_row)

def get_bishop_moves(pos, occupied_squares):
    """Generates all legal bishop moves."""
    col, row = pos
    # Directions: up-right, up-left, down-right, down-left
    for dc, dr in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
        for i in range(1, 8):
            new_pos = (col + i * dc, row + i * dr)
            if not (0 <= new_pos[0] < 8 and 0 <= new_pos[1] < 8):
                break
            yield new_pos
            if new_pos in occupied_squares:
                break

memo = {}

def find_forced_mate(wk_pos, wb_pos, bk_pos, turn, depth):
    """
    Minimax search function to find a forced mate.
    'depth' here refers to plies (half-moves).
    Returns a winning move if one is found, otherwise None.
    """
    state_key = (wk_pos, wb_pos, bk_pos, turn, depth)
    if state_key in memo:
        return memo[state_key]

    # Get squares attacked by white
    attacked_by_white = set(get_king_moves(wk_pos))
    for move in get_bishop_moves(wb_pos, {wk_pos, bk_pos}):
        attacked_by_white.add(move)

    # Black's turn logic
    if turn == 'B':
        # Get black king's legal moves
        bk_legal_moves = []
        for move in get_king_moves(bk_pos):
            if is_white_square(move) and move not in attacked_by_white:
                bk_legal_moves.append(move)

        is_in_check = bk_pos in attacked_by_white

        if not bk_legal_moves:
            # If no legal moves, it's either checkmate or stalemate
            memo[state_key] = 'win' if is_in_check else 'draw'
            return memo[state_key]
        
        if depth == 0:
            memo[state_key] = None
            return None

        # Assume opponent can escape unless all replies lead to a loss
        for bk_move in bk_legal_moves:
            # Recursively check if white can force a mate after this black move
            if find_forced_mate(wk_pos, wb_pos, bk_move, 'W', depth - 1) is None:
                memo[state_key] = None  # Black found an escape
                return None
        
        # All of black's replies lead to a position where white can win
        memo[state_key] = 'win'
        return 'win'

    # White's turn logic
    if turn == 'W':
        if depth == 0:
            memo[state_key] = None
            return None

        # Try all white king moves
        bk_danger_squares = set(get_king_moves(bk_pos))
        for wk_move in get_king_moves(wk_pos):
            if wk_move not in bk_danger_squares and wk_move != wb_pos:
                if find_forced_mate(wk_move, wb_pos, bk_pos, 'B', depth - 1) == 'win':
                    memo[state_key] = ('K', wk_move)
                    return ('K', wk_move)

        # Try all white bishop moves
        for wb_move in get_bishop_moves(wb_pos, {wk_pos, bk_pos}):
            if find_forced_mate(wk_pos, wb_move, bk_pos, 'B', depth - 1) == 'win':
                memo[state_key] = ('B', wb_move)
                return ('B', wb_move)
        
        memo[state_key] = None
        return None


def solve():
    """
    Iteratively searches for the shortest mate.
    """
    wk_initial = alg_to_coord('d2')
    wb_initial = alg_to_coord('e2')
    bk_initial = alg_to_coord('d5')

    # Max search depth for mate in N moves is 2*N-1 plies
    for d in range(1, 15): # Look for mate in 1 up to mate in 7
        # Clear memo for new depth search, but could be optimized
        global memo
        memo.clear()
        
        # Search for a mate in `d` moves (which is `2d-1` plies from White's first move)
        # We start with White's turn, so total plies = 2d - 1
        result = find_forced_mate(wk_initial, wb_initial, bk_initial, 'W', 2 * d - 1)
        
        if result:
            piece, move_coord = result
            move_str = piece_to_alg(piece, move_coord)
            print(f"{move_str}, {d}")
            return

if __name__ == '__main__':
    solve()