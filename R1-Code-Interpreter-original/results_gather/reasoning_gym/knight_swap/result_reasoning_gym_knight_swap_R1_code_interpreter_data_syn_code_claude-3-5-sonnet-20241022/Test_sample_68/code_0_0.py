from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 4 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_position(piece, board):
    positions = []
    for i in range(4):
        for j in range(4):
            if board[i][j] == piece:
                positions.append((i, j))
    return positions

def is_target_state(board):
    white_positions = {(0,2), (2,2)}  # Original white positions
    black_positions = {(1,3), (2,0)}  # Original black positions
    
    current_whites = set(get_position('w', board))
    current_blacks = set(get_position('B', board))
    
    return current_whites == black_positions and current_blacks == white_positions

def convert_move(from_pos, to_pos):
    cols = 'ABCD'
    rows = '123'
    from_notation = f"{cols[from_pos[1]]}{rows[from_pos[0]]}"
    to_notation = f"{cols[to_pos[1]]}{rows[to_pos[0]]}"
    return from_notation, to_notation

def solve_knight_swap():
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