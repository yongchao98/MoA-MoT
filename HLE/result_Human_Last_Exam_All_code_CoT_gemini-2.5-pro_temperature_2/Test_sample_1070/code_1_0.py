import collections

def to_coords(alg):
    """Converts algebraic notation like 'd5' to (x, y) coordinates (0-indexed)."""
    col = ord(alg[0]) - ord('a')
    row = int(alg[1]) - 1
    return col, row

def to_alg(coords):
    """Converts (x, y) coordinates to algebraic notation."""
    col_char = chr(coords[0] + ord('a'))
    row_char = str(coords[1] + 1)
    return col_char + row_char

def is_on_board(coords):
    """Checks if coordinates are on the 8x8 board."""
    x, y = coords
    return 0 <= x < 8 and 0 <= y < 8

def is_white_square_coords(coords):
    """Checks if a square is white. (x,y) is white if x+y is odd."""
    x, y = coords
    return (x + y) % 2 != 0

def get_piece_alg(piece_char, target_coords):
    """Formats a move in algebraic notation, e.g., 'Bc4'."""
    return piece_char + to_alg(target_coords)

def get_king_moves(coords):
    """Generates all 8 potential king move destinations."""
    x, y = coords
    moves = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            new_coords = (x + dx, y + dy)
            if is_on_board(new_coords):
                moves.append(new_coords)
    return moves

def get_bishop_moves(bishop_coords, occupied_squares):
    """Generates all legal bishop moves."""
    moves = []
    x, y = bishop_coords
    directions = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    for dx, dy in directions:
        for i in range(1, 8):
            new_coords = (x + i * dx, y + i * dy)
            if not is_on_board(new_coords):
                break
            if new_coords in occupied_squares:
                break
            moves.append(new_coords)
    return moves

def is_square_attacked_by_white(square_coords, wk_coords, wb_coords):
    """Checks if a square is attacked by the white king or bishop."""
    # King attack
    if abs(square_coords[0] - wk_coords[0]) <= 1 and abs(square_coords[1] - wk_coords[1]) <= 1:
        return True
    # Bishop attack
    if abs(square_coords[0] - wb_coords[0]) == abs(square_coords[1] - wb_coords[1]):
        # Check for blocking pieces (only the white king could block)
        dx = 1 if square_coords[0] > wb_coords[0] else -1
        dy = 1 if square_coords[1] > wb_coords[1] else -1
        cx, cy = wb_coords[0] + dx, wb_coords[1] + dy
        while (cx, cy) != square_coords:
            if (cx, cy) == wk_coords:
                return False # Path is blocked by the white king
            cx += dx
            cy += dy
        return True
    return False

def get_black_king_legal_moves(bk_coords, wk_coords, wb_coords):
    """Gets legal moves for the black king (must be to a white, unattacked square)."""
    legal_moves = []
    potential_moves = get_king_moves(bk_coords)
    for move in potential_moves:
        # King cannot move into adjacency with the other king
        if abs(move[0] - wk_coords[0]) <= 1 and abs(move[1] - wk_coords[1]) <= 1:
            continue
        if is_white_square_coords(move) and not is_square_attacked_by_white(move, wk_coords, wb_coords):
            legal_moves.append(move)
    return legal_moves

def solve():
    """Finds the shortest checkmate using BFS."""
    wk_start = to_coords('d2')
    wb_start = to_coords('e2')
    bk_start = to_coords('d5')
    
    initial_state = (wk_start, wb_start, bk_start)
    
    # queue stores: (state, path_of_moves)
    queue = collections.deque([(initial_state, [])])
    visited = {initial_state}

    while queue:
        (wk_pos, wb_pos, bk_pos), path = queue.popleft()

        # It's White's turn
        
        # 1. Generate White King moves
        wk_moves = get_king_moves(wk_pos)
        for new_wk_pos in wk_moves:
             # King can't move adjacent to enemy king
            if abs(new_wk_pos[0] - bk_pos[0]) <= 1 and abs(new_wk_pos[1] - bk_pos[1]) <= 1:
                continue
            # King cannot move to where the bishop is
            if new_wk_pos == wb_pos:
                continue
            
            # Form move and new state
            move_alg = get_piece_alg('K', new_wk_pos)
            current_move_path = path + [move_alg]

            # Check for mate
            is_check = is_square_attacked_by_white(bk_pos, new_wk_pos, wb_pos)
            bk_legal_moves = get_black_king_legal_moves(bk_pos, new_wk_pos, wb_pos)

            if is_check and not bk_legal_moves:
                num_moves = len(current_move_path) // 2 + 1
                first_move = current_move_path[0]
                print(f"{first_move}, {num_moves}")
                return

            # If not mate, queue black's responses
            for new_bk_pos in bk_legal_moves:
                new_state = (new_wk_pos, wb_pos, new_bk_pos)
                if new_state not in visited:
                    bk_move_alg = get_piece_alg('K', new_bk_pos)
                    visited.add(new_state)
                    queue.append((new_state, current_move_path + [bk_move_alg]))

        # 2. Generate White Bishop moves
        wb_moves = get_bishop_moves(wb_pos, {wk_pos, bk_pos})
        for new_wb_pos in wb_moves:
            # Form move and new state
            move_alg = get_piece_alg('B', new_wb_pos)
            current_move_path = path + [move_alg]
            
            # Check for mate
            is_check = is_square_attacked_by_white(bk_pos, wk_pos, new_wb_pos)
            bk_legal_moves = get_black_king_legal_moves(bk_pos, wk_pos, new_wb_pos)
            
            if is_check and not bk_legal_moves:
                num_moves = len(current_move_path) // 2 + 1
                first_move = current_move_path[0]
                print(f"{first_move}, {num_moves}")
                return
            
            # If not mate, queue black's responses
            for new_bk_pos in bk_legal_moves:
                new_state = (wk_pos, new_wb_pos, new_bk_pos)
                if new_state not in visited:
                    bk_move_alg = get_piece_alg('K', new_bk_pos)
                    visited.add(new_state)
                    queue.append((new_state, current_move_path + [bk_move_alg]))
    
    print("No checkmate found.")

solve()
<<<Bc4, 6>>>