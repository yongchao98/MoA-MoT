from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [(row-2,col+1), (row-2,col-1), (row+2,col+1), (row+2,col-1),
             (row-1,col+2), (row-1,col-2), (row+1,col+2), (row+1,col-2)]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]  # Fixed bounds

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_position_name(pos):
    row, col = pos
    return chr(ord('A') + col) + str(3 - row)  # Adjusted for chess notation

def solve_knight_swap():
    # Initial board setup (0-based indexing, top to bottom)
    initial_board = [
        ['.', 'B', '.', ' '],  # row 3 (index 0)
        [' ', ' ', '.', 'w'],  # row 2 (index 1)
        ['.', 'w', ' ', 'B']   # row 1 (index 2)
    ]
    
    # Target positions for each color (row, col)
    target_positions = {
        'w': {(0,1), (2,3)},  # B3, D1 positions
        'B': {(2,1), (1,3)}   # B1, D2 positions
    }
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        # Check if we reached target state
        current_positions = {'w': set(), 'B': set()}
        for i in range(3):
            for j in range(4):
                if board[i][j] in ['w', 'B']:
                    current_positions[board[i][j]].add((i,j))
        
        if (current_positions['w'] == target_positions['w'] and 
            current_positions['B'] == target_positions['B']):
            return moves
        
        # Generate next moves
        piece = 'w' if is_white_turn else 'B'
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    for new_i, new_j in get_knight_moves((i,j)):
                        if board[new_i][new_j] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[i][j] = '.'
                            new_board[new_i][new_j] = piece
                            move = f"{piece},{get_position_name((i,j))},{get_position_name((new_i,new_j))}"
                            queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")