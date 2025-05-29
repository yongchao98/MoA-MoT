from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def convert_to_chess_notation(x, y):
    # Convert from 0-based array indices to chess notation
    return f"{chr(65+x)}{3-y}"

def is_final_state(board):
    # White knights should be at C1, D1 (board[2][2] and board[2][3])
    # Black knights should be at B3, D2 (board[0][1] and board[1][3])
    return (board[2][2] == 'w' and board[2][3] == 'w' and 
            board[0][1] == 'B' and board[1][3] == 'B')

def find_solution():
    # Initial board state (0-based indexing)
    # row 3 (index 0), row 2 (index 1), row 1 (index 2)
    initial_board = [
        ['.', 'w', '.', '.'],  # row 3
        ['.', '.', '.', 'w'],  # row 2
        ['.', '.', 'B', 'B']   # row 1
    ]
    
    # Queue stores: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, white_turn = queue.popleft()
        
        # Limit search depth
        if len(moves) > 8:
            continue
        
        # Find all pieces for current player
        current_pieces = []
        for y in range(3):
            for x in range(4):
                if current_board[y][x] == ('w' if white_turn else 'B'):
                    current_pieces.append((x, y))
        
        # Try all possible moves for each piece
        for piece_pos in current_pieces:
            for new_pos in get_knight_moves(piece_pos):
                # Check if target square is empty
                if current_board[new_pos[1]][new_pos[0]] == '.':
                    # Create new board with the move
                    new_board = copy.deepcopy(current_board)
                    new_board[piece_pos[1]][piece_pos[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = 'w' if white_turn else 'B'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        
                        # Record move in chess notation
                        from_square = convert_to_chess_notation(piece_pos[0], piece_pos[1])
                        to_square = convert_to_chess_notation(new_pos[0], new_pos[1])
                        move = f"{'w' if white_turn else 'B'},{from_square},{to_square}"
                        new_moves = moves + [move]
                        
                        # Check if we reached the target state
                        if is_final_state(new_board):
                            return new_moves
                        
                        # Add new state to queue
                        queue.append((new_board, new_moves, not white_turn))
    
    return "No"

# Run the solution
result = find_solution()
print(result)