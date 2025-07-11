import collections

def is_light(square_str):
    """Checks if a square is light (white). a1 is dark."""
    if not square_str or len(square_str) != 2:
        return False
    col = ord(square_str[0]) - ord('a')
    row = int(square_str[1]) - 1
    return (col + row) % 2 != 0

def get_king_moves(pos):
    """Get all king moves from a position."""
    col, row = pos
    moves = set()
    for dc in [-1, 0, 1]:
        for dr in [-1, 0, 1]:
            if dc == 0 and dr == 0:
                continue
            nc, nr = col + dc, row + dr
            if 0 <= nc < 8 and 0 <= nr < 8:
                moves.add((nc, nr))
    return moves

def to_alg(pos):
    """Convert (col, row) to algebraic notation."""
    col, row = pos
    return chr(ord('a') + col) + str(row + 1)

def to_pos(alg):
    """Convert algebraic notation to (col, row)."""
    return (ord(alg[0]) - ord('a'), int(alg[1]) - 1)

def get_bishop_attacks(pos, board_state):
    """Get all squares a bishop at pos attacks."""
    col, row = pos
    attacks = set()
    # 4 directions
    for dc in [-1, 1]:
        for dr in [-1, 1]:
            for i in range(1, 8):
                nc, nr = col + i * dc, row + i * dr
                if 0 <= nc < 8 and 0 <= nr < 8:
                    attacks.add((nc, nr))
                    if (nc,nr) in board_state: # Stop if piece is hit
                        break
                else:
                    break
    return attacks


# We will solve this by a Breadth-First Search (BFS) to find the shortest path to mate.
# A state is (wk_pos, wb_pos, bk_pos, turn), where turn is 'w' or 'b'
initial_state = (to_pos('d2'), to_pos('e2'), to_pos('d5'), 'w')

# queue stores (state, path_list)
q = collections.deque([(initial_state, [])])
# visited prevents cycles and redundant computations
visited = {initial_state}

mate_path = None

while q:
    current_state, path = q.popleft()
    wk_pos, wb_pos, bk_pos, turn = current_state

    # For white's turn, we generate all possible moves for white
    if turn == 'w':
        board_state = {wk_pos: 'WK', wb_pos: 'WB', bk_pos: 'BK'}
        
        # King moves
        for next_wk_pos in get_king_moves(wk_pos):
            # King can't move to a square attacked by the other king
            if next_wk_pos not in get_king_moves(bk_pos) and next_wk_pos != wb_pos:
                move_str = "K" + to_alg(next_wk_pos)
                new_state = (next_wk_pos, wb_pos, bk_pos, 'b')
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((new_state, path + [move_str]))
        
        # Bishop moves
        for next_wb_pos in get_bishop_attacks(wb_pos, board_state):
             if next_wb_pos != wk_pos:
                move_str = "B" + to_alg(next_wb_pos)
                new_state = (wk_pos, next_wb_pos, bk_pos, 'b')
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((new_state, path + [move_str]))

    # For black's turn, we find legal moves and check for mate/stalemate
    elif turn == 'b':
        # What squares does white control?
        white_king_attacks = get_king_moves(wk_pos)
        white_bishop_attacks = get_bishop_attacks(wb_pos, {wk_pos: 'WK', bk_pos: 'BK'})
        white_controlled_squares = white_king_attacks.union(white_bishop_attacks)

        # Find legal moves for black king
        legal_bk_moves = []
        for next_bk_pos in get_king_moves(bk_pos):
            # 'Cursed' to only move on light squares and can't move into check
            if is_light(to_alg(next_bk_pos)) and next_bk_pos not in white_controlled_squares:
                 legal_bk_moves.append(next_bk_pos)
        
        is_in_check = bk_pos in white_controlled_squares

        if not legal_bk_moves:
            if is_in_check:
                # CHECKMATE FOUND
                mate_path = path
                break 
            else:
                # STALEMATE (draw), ignore this path
                continue
        
        # Add next possible black moves to queue
        for next_bk_pos in legal_bk_moves:
             move_str = "K" + to_alg(next_bk_pos)
             new_state = (wk_pos, wb_pos, next_bk_pos, 'w')
             # Visited check here might be complex because black could choose non-optimal moves
             # but for BFS shortest path, we assume optimal play. So we only need to visit a state once.
             if new_state not in visited:
                visited.add(new_state)
                q.append((new_state, path + [move_str]))


if mate_path:
    first_move = mate_path[0]
    num_moves = (len(mate_path) + 1) // 2
    print(f"{first_move}, {num_moves}")
else:
    print("Could not find a checkmate sequence.")
