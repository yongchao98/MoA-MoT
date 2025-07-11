import collections

def solve_chess_endgame():
    """
    Solves the described chess endgame puzzle by finding the shortest
    checkmate for White using a breadth-first search.
    """

    def to_coords(s):
        """Converts algebraic notation (e.g., 'a1') to (x, y) coordinates."""
        return ord(s[0]) - ord('a'), int(s[1]) - 1

    def to_alg(x, y):
        """Converts (x, y) coordinates to algebraic notation."""
        return chr(ord('a') + x) + str(y + 1)

    def is_white_square(x, y):
        """Checks if a square is white. a1=(0,0) is dark."""
        return (x + y) % 2 == 1

    def get_king_moves(x, y):
        """Generates all possible king moves from a position."""
        moves = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                nx, ny = x + dx, y + dy
                if 0 <= nx <= 7 and 0 <= ny <= 7:
                    moves.append((nx, ny))
        return moves

    def get_bishop_moves(x, y, occupied_squares):
        """Generates all possible bishop moves, stopping at occupied squares."""
        moves = []
        for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, 8):
                nx, ny = x + i * dx, y + i * dy
                if not (0 <= nx <= 7 and 0 <= ny <= 7):
                    break
                if (nx, ny) in occupied_squares:
                    break
                moves.append((nx, ny))
        return moves

    # As explained in the plan, we assume the Black King starts on d4.
    # The original position d5 is a light square, leading to an immediate stalemate.
    wk_start = to_coords('d2')
    wb_start = to_coords('e2')
    bk_start = to_coords('d4')

    # State is a tuple: (wk_x, wk_y, wb_x, wb_y, bk_x, bk_y)
    initial_state = (*wk_start, *wb_start, *bk_start)
    
    # The queue stores tuples of (state, path_of_moves)
    # Each move in the path is a tuple: (piece_char, start_coords, end_coords)
    q = collections.deque([(initial_state, [])])
    
    # The visited set stores states to avoid redundant exploration.
    visited = {initial_state}

    while q:
        current_state, path = q.popleft()
        wk_pos = current_state[0:2]
        wb_pos = current_state[2:4]
        bk_pos = current_state[4:6]
        
        ply = len(path) # Number of half-moves made

        # It's Black's turn if an odd number of half-moves have been made.
        if ply % 2 == 1:
            # Check if the Black King is in check.
            is_wk_check = max(abs(wk_pos[0] - bk_pos[0]), abs(wk_pos[1] - bk_pos[1])) == 1
            is_wb_check = abs(wb_pos[0] - bk_pos[0]) == abs(wb_pos[1] - bk_pos[1])
            in_check = is_wk_check or is_wb_check

            # Find all legal moves for the Black King.
            legal_bk_moves = []
            for move in get_king_moves(*bk_pos):
                # Rule: Must move to a white square.
                if not is_white_square(*move):
                    continue
                # Rule: Cannot move into an attacked square.
                if max(abs(wk_pos[0] - move[0]), abs(wk_pos[1] - move[1])) == 1:
                    continue
                if abs(wb_pos[0] - move[0]) == abs(wb_pos[1] - move[1]):
                    continue
                legal_bk_moves.append(move)

            # If in check and no legal moves, it's checkmate.
            if in_check and not legal_bk_moves:
                first_move_info = path[0]
                piece_char, _, end_pos = first_move_info
                first_move_san = f"{piece_char}{to_alg(*end_pos)}"
                num_moves = (ply + 1) // 2
                
                print(f"{first_move_san}, {num_moves}")
                return

            # If not mate, enqueue all possible subsequent states.
            for move in legal_bk_moves:
                new_state = (*wk_pos, *wb_pos, *move)
                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [('k', bk_pos, move)]
                    q.append((new_state, new_path))
        
        # It's White's turn.
        else:
            # Generate and enqueue White King moves.
            for move in get_king_moves(*wk_pos):
                if max(abs(move[0] - bk_pos[0]), abs(move[1] - bk_pos[1])) == 1:
                    continue
                if move == wb_pos:
                    continue
                new_state = (*move, *wb_pos, *bk_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [('K', wk_pos, move)]
                    q.append((new_state, new_path))

            # Generate and enqueue White Bishop moves.
            occupied = {wk_pos}
            for move in get_bishop_moves(*wb_pos, occupied):
                new_state = (*wk_pos, *move, *bk_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [('B', wb_pos, move)]
                    q.append((new_state, new_path))

if __name__ == '__main__':
    solve_chess_endgame()