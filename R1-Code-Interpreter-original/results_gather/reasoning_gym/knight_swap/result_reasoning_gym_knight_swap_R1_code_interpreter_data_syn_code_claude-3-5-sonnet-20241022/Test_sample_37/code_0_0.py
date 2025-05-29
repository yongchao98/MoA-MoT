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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are where black knights started and vice versa
    initial_whites = {(2,0), (2,3)}
    initial_blacks = {(2,1), (0,0)}
    
    current_whites = set()
    current_blacks = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                current_whites.add((i,j))
            elif board[i][j] == 'B':
                current_blacks.add((i,j))
    
    return current_whites == initial_blacks and current_blacks == initial_whites

def convert_move(from_pos, to_pos):
    cols = 'ABCD'
    rows = '123'
    from_notation = cols[from_pos[1]] + rows[2-from_pos[0]]
    to_notation = cols[to_pos[1]] + rows[2-to_pos[0]]
    return from_notation + "," + to_notation

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['w', 'B', '.', 'w'],
        ['.', '.', '.', '.'],
        ['B', '.', '.', '.']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(board):
            return moves
        
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('w' if is_white_turn else 'B'):
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            valid_moves = get_knight_moves(piece_pos)
            for new_pos in valid_moves:
                if board[new_pos[0]][new_pos[1]] == '.':  # If target is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[new_pos[0]][new_pos[1]] = new_board[piece_pos[0]][piece_pos[1]]
                    new_board[piece_pos[0]][piece_pos[1]] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        piece_color = 'w' if is_white_turn else 'B'
                        move = piece_color + "," + convert_move(piece_pos, new_pos)
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

# Solve and print result
result = solve_knight_swap()
if result is None:
    print("<<<No>>>")
else:
    print(f"<<<{result}>>>")