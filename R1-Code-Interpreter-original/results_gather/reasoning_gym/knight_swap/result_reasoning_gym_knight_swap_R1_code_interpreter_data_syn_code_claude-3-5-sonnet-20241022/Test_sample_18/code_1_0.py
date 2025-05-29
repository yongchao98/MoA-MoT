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
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at B3 and B2 positions
    # and black knights are at A1 and B1 positions
    return (board[0][1] == 'w' and board[1][1] == 'w' and 
            board[2][0] == 'B' and board[2][1] == 'B')

def get_valid_moves(board, is_white_turn):
    piece = 'w' if is_white_turn else 'B'
    moves = []
    
    # Find all pieces of current player
    for i in range(3):
        for j in range(4):
            if board[i][j] == piece:
                start = (i, j)
                # Get all possible knight moves
                for end in get_knight_moves(start):
                    if board[end[0]][end[1]] == '.':  # If target square is empty
                        moves.append((start, end))
    return moves

def make_move(board, start, end):
    new_board = [row[:] for row in board]
    piece = new_board[start[0]][start[1]]
    new_board[start[0]][start[1]] = '.'
    new_board[end[0]][end[1]] = piece
    return new_board

def solve_knight_swap():
    # Initial board state (3x4)
    initial_board = [
        ['.', 'B', '.', '.'],
        ['B', '.', '.', '.'],
        ['w', 'w', '.', '.']
    ]
    
    queue = deque([(initial_board, True, [])])
    seen = set()
    
    while queue:
        board, is_white_turn, moves = queue.popleft()
        board_str = board_to_string(board) + str(is_white_turn)
        
        if board_str in seen:
            continue
        seen.add(board_str)
        
        if is_target_reached(board):
            return moves
        
        valid_moves = get_valid_moves(board, is_white_turn)
        for start, end in valid_moves:
            new_board = make_move(board, start, end)
            # Convert coordinates to chess notation (A1 is bottom-left)
            move_notation = f"{'w' if is_white_turn else 'B'},{chr(65+start[1])}{3-start[0]},{chr(65+end[1])}{3-end[0]}"
            new_moves = moves + [move_notation]
            queue.append((new_board, not is_white_turn, new_moves))
    
    return None

# Solve and print result
result = solve_knight_swap()
if result is None:
    print("<<<No>>>")
else:
    print(f"<<<{result}>>>")