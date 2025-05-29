from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]  # Fixed board dimensions

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_position(piece, board):
    positions = []
    for i in range(3):  # Fixed board dimensions
        for j in range(4):
            if board[i][j] == piece:
                positions.append((i, j))
    return positions

def is_target_state(board):
    # Check if white knights are at original black positions
    white_pos = set(get_position('w', board))
    black_pos = set(get_position('B', board))
    original_white = {(0,2), (2,2)}  # C3, C1
    original_black = {(1,3), (2,0)}  # D2, A1
    return white_pos == original_black and black_pos == original_white

def convert_move(from_pos, to_pos):
    cols = 'ABCD'
    rows = '321'  # Reversed to match the board representation
    from_notation = f"{cols[from_pos[1]]}{rows[from_pos[0]]}"
    to_notation = f"{cols[to_pos[1]]}{rows[to_pos[0]]}"
    return from_notation, to_notation

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', '.', 'w', '.'],
        ['.', '.', '.', 'B'],
        ['B', '.', 'w', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(current_board):
            return moves
        
        piece = 'w' if is_white_turn else 'B'
        positions = get_position(piece, current_board)
        
        for pos in positions:
            valid_moves = get_knight_moves(pos)
            for new_pos in valid_moves:
                if current_board[new_pos[0]][new_pos[1]] == '.':
                    new_board = [row[:] for row in current_board]
                    new_board[new_pos[0]][new_pos[1]] = piece
                    new_board[pos[0]][pos[1]] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        from_not, to_not = convert_move(pos, new_pos)
                        new_moves = moves + [f"{piece},{from_not},{to_not}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")