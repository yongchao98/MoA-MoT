from collections import deque

def get_knight_moves(pos, size=4):
    moves = [
        (2,1), (2,-1), (-2,1), (-2,-1),
        (1,2), (1,-2), (-1,2), (-1,-2)
    ]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < size and 0 <= new_y < size:
            valid.append((new_x, new_y))
    return valid

def pos_to_coord(pos):
    return f"{chr(65+pos[1])}{4-pos[0]}"  # Adjusted for correct coordinate system

def coord_to_pos(coord):
    return (4-int(coord[1]), ord(coord[0])-65)

def get_pieces_positions(board):
    positions = {'w': [], 'B': []}
    for i in range(4):
        for j in range(4):
            if board[i][j] in 'wB':
                positions[board[i][j]].append((i, j))
    return positions

def is_target_reached(positions):
    # White knights should be at original Black positions and vice versa
    w_target = {(0,2), (1,3)}  # Original Black positions
    b_target = {(1,2), (0,3)}  # Original White positions
    return set(positions['w']) == w_target and set(positions['B']) == b_target

def make_move(board, from_pos, to_pos, piece):
    new_board = [row[:] for row in board]
    new_board[from_pos[0]][from_pos[1]] = '.'
    new_board[to_pos[0]][to_pos[1]] = piece
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_knight_swap():
    # Initialize board (0,0 is top-left)
    initial_board = [
        ['.', '.', 'B', 'w'],
        ['.', '.', 'w', 'B'],
        ['.', '.', '.', '.'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        positions = get_pieces_positions(board)
        
        if is_target_reached(positions):
            return moves
        
        current_piece = 'w' if is_white_turn else 'B'
        for piece_pos in positions[current_piece]:
            for new_pos in get_knight_moves(piece_pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    new_board = make_move(board, piece_pos, new_pos, current_piece)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_move = f"{current_piece},{pos_to_coord(piece_pos)},{pos_to_coord(new_pos)}"
                        queue.append((new_board, moves + [new_move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")