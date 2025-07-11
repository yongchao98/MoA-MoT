import collections

def solve():
    """
    Solves the cursed king chess puzzle using Breadth-First Search.
    """
    COLS = 'abcdefgh'
    ROWS = '12345678'

    def to_algebraic(pos):
        col, row = pos
        return COLS[col] + ROWS[row]

    def is_on_board(pos):
        col, row = pos
        return 0 <= col < 8 and 0 <= row < 8

    def is_white_square(pos):
        col, row = pos
        return (col + row) % 2 != 0

    # --- Move Generation and Game Logic ---

    def get_bishop_attacks(wb_pos, occupied_squares):
        attacks = set()
        for dc, dr in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, 8):
                pos = (wb_pos[0] + i * dc, wb_pos[1] + i * dr)
                if not is_on_board(pos):
                    break
                attacks.add(pos)
                if pos in occupied_squares:
                    break
        return attacks

    def get_king_attacks(k_pos):
        attacks = set()
        for dc in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if dc == 0 and dr == 0:
                    continue
                pos = (k_pos[0] + dc, k_pos[1] + dr)
                if is_on_board(pos):
                    attacks.add(pos)
        return attacks
        
    def get_white_moves(wk_pos, bk_pos, wb_pos):
        moves = []
        # King moves
        for dc in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if dc == 0 and dr == 0:
                    continue
                next_pos = (wk_pos[0] + dc, wk_pos[1] + dr)
                if is_on_board(next_pos) and max(abs(next_pos[0] - bk_pos[0]), abs(next_pos[1] - bk_pos[1])) > 1:
                     moves.append(('K', next_pos))
    
        # Bishop moves
        for dc, dr in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, 8):
                next_pos = (wb_pos[0] + i * dc, wb_pos[1] + i * dr)
                if not is_on_board(next_pos):
                    break
                if next_pos == wk_pos or next_pos == bk_pos:
                    break
                moves.append(('B', next_pos))
        return moves

    def get_black_king_moves(bk_pos, wk_pos, wb_pos):
        moves = []
        white_king_attacks = get_king_attacks(wk_pos)
        # Note: black king cannot block bishop check on itself
        white_bishop_attacks = get_bishop_attacks(wb_pos, {wk_pos}) 
        attacked_by_white = white_king_attacks | white_bishop_attacks

        for dc in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if dc == 0 and dr == 0:
                    continue
                next_pos = (bk_pos[0] + dc, bk_pos[1] + dr)
                if is_on_board(next_pos) and is_white_square(next_pos) and next_pos not in attacked_by_white:
                    moves.append(next_pos)
        return moves
        
    def is_checkmate(wk_pos, bk_pos, wb_pos):
        # Is black king in check from bishop?
        # A king cannot deliver check since kings cannot be adjacent
        if bk_pos not in get_bishop_attacks(wb_pos, {wk_pos}):
            return False
        # Does black have any legal moves?
        if len(get_black_king_moves(bk_pos, wk_pos, wb_pos)) == 0:
            return True
        return False

    # --- BFS Solver ---

    # Initial positions: BK d5, WK d2, WB e2
    initial_state = ((3, 1), (3, 4), (4, 1), 'w') # (wk_pos, bk_pos, wb_pos, turn)

    # Queue stores tuples of (state, path_list)
    queue = collections.deque([(initial_state, [])]) 
    # Visited set stores states to prevent cycles
    visited = {initial_state}

    while queue:
        (wk_pos, bk_pos, wb_pos, turn), path = queue.popleft()

        if turn == 'w':
            # Generate all legal moves for White
            white_moves = get_white_moves(wk_pos, bk_pos, wb_pos)
            for piece, next_piece_pos in white_moves:
                move_notation = piece + to_algebraic(next_piece_pos)
                new_path = path + [move_notation]

                if piece == 'K':
                    next_wk_pos, next_wb_pos = next_piece_pos, wb_pos
                else: # Bishop
                    next_wk_pos, next_wb_pos = wk_pos, next_piece_pos
                
                # Check for mate after White's move
                if is_checkmate(next_wk_pos, bk_pos, next_wb_pos):
                    first_move = new_path[0]
                    num_white_moves = (len(new_path) + 1) // 2
                    print(f"{first_move}, {num_white_moves}")
                    # Also print the final checkmating sequence
                    final_equation = []
                    for i, move in enumerate(new_path):
                        if i % 2 == 0:
                            final_equation.append(f"{i//2 + 1}. {move}")
                        else:
                            final_equation[-1] += f" {move}"
                    print("Full sequence:", ' '.join(final_equation))
                    return

                next_state = (next_wk_pos, bk_pos, next_wb_pos, 'b')
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append((next_state, new_path))
        
        else: # turn == 'b'
            # Generate legal moves for Black King
            black_moves = get_black_king_moves(bk_pos, wk_pos, wb_pos)
            for next_bk_pos in black_moves:
                move_notation = 'k' + to_algebraic(next_bk_pos) # unofficial notation
                new_path = path + [move_notation]
                next_state = (wk_pos, next_bk_pos, wb_pos, 'w')
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append((next_state, new_path))

solve()