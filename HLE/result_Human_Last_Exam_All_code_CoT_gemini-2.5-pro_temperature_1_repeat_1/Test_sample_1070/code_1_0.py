import collections

def solve_chess_endgame():
    """
    Solves the specified chess endgame with a cursed black king using Breadth-First Search.
    Finds the shortest checkmate for White.
    """

    # --- Helper Functions ---
    COLS = "abcdefgh"
    ROWS = "12345678"

    def to_coords(notation):
        col = COLS.index(notation[0])
        row = ROWS.index(notation[1])
        return col, row

    def to_alg(col, row):
        return f"{COLS[col]}{ROWS[row]}"

    def is_valid(col, row):
        return 0 <= col < 8 and 0 <= row < 8

    def is_white_square(col, row):
        return (col + row) % 2 != 0

    def get_move_notation(prev_state, curr_state):
        p_wkc, p_wkr, p_wbc, p_wbr, _, _ = prev_state
        c_wkc, c_wkr, c_wbc, c_wbr, _, _ = curr_state

        if (p_wkc, p_wkr) != (c_wkc, c_wkr):
            return f"K{to_alg(c_wkc, c_wkr)}"
        if (p_wbc, p_wbr) != (c_wbc, c_wbr):
            return f"B{to_alg(c_wbc, c_wbr)}"
        return ""

    def is_attacked(sq_col, sq_row, wk_c, wk_r, wb_c, wb_r, blocker_c, blocker_r):
        # Attacked by White King
        if abs(sq_col - wk_c) <= 1 and abs(sq_row - wk_r) <= 1:
            return True
        
        # Attacked by White Bishop
        if abs(sq_col - wb_c) == abs(sq_row - wb_r):
            dc = 1 if sq_col > wb_c else -1
            dr = 1 if sq_row > wb_r else -1
            cc, cr = wb_c + dc, wb_r + dr
            is_blocked = False
            while (cc, cr) != (sq_col, sq_row):
                if (cc, cr) == (wk_c, wk_r) or (cc, cr) == (blocker_c, blocker_r):
                    is_blocked = True
                    break
                cc += dc
                cr += dr
            if not is_blocked:
                return True
        return False

    # --- BFS Initialization ---
    initial_state = (*to_coords('d2'), *to_coords('e2'), *to_coords('d5'))
    queue = collections.deque([(initial_state, [initial_state])])
    visited = {initial_state}

    # --- Main BFS Loop ---
    while queue:
        current_state, path = queue.popleft()
        wk_c, wk_r, wb_c, wb_r, bk_c, bk_r = current_state
        is_white_turn = (len(path) - 1) % 2 == 0

        if is_white_turn:
            # Generate White King moves
            for dc, dr in [(0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1), (-1,0), (-1,1)]:
                nwk_c, nwk_r = wk_c + dc, wk_r + dr
                if is_valid(nwk_c, nwk_r) and (abs(nwk_c - bk_c) > 1 or abs(nwk_r - bk_r) > 1):
                    next_state = (nwk_c, nwk_r, wb_c, wb_r, bk_c, bk_r)
                    if next_state not in visited:
                        visited.add(next_state)
                        queue.append((next_state, path + [next_state]))

            # Generate White Bishop moves
            for dc, dr in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                for i in range(1, 8):
                    nwb_c, nwb_r = wb_c + i * dc, wb_r + i * dr
                    if not is_valid(nwb_c, nwb_r): break
                    
                    next_state = (wk_c, wk_r, nwb_c, nwb_r, bk_c, bk_r)
                    if next_state not in visited:
                        visited.add(next_state)
                        queue.append((next_state, path + [next_state]))
                    
                    if (nwb_c, nwb_r) == (wk_c, wk_r) or (nwb_c, nwb_r) == (bk_c, bk_r):
                        break
        else:  # Black's turn
            is_in_check = is_attacked(bk_c, bk_r, wk_c, wk_r, wb_c, wb_r, bk_c, bk_r)
            legal_black_moves = []
            for dc, dr in [(0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1), (-1,0), (-1,1)]:
                nbk_c, nbk_r = bk_c + dc, bk_r + dr
                if (is_valid(nbk_c, nbk_r) and is_white_square(nbk_c, nbk_r) and
                    (abs(nbk_c - wk_c) > 1 or abs(nbk_r - wk_r) > 1) and
                    not is_attacked(nbk_c, nbk_r, wk_c, wk_r, wb_c, wb_r, nbk_c, nbk_r)):
                    legal_black_moves.append((nbk_c, nbk_r))

            if is_in_check and not legal_black_moves:
                first_move = get_move_notation(path[0], path[1])
                num_moves = (len(path) - 1) // 2 + 1
                print(f"{first_move}, {num_moves}")
                return

            for nbk_c, nbk_r in legal_black_moves:
                next_state = (wk_c, wk_r, wb_c, wb_r, nbk_c, nbk_r)
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append((next_state, path + [next_state]))

if __name__ == '__main__':
    solve_chess_endgame()