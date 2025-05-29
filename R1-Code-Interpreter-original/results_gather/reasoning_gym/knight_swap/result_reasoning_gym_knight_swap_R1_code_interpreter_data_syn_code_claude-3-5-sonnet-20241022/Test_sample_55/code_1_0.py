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
    # Adjust for 0-based indexing
    return [(r-1,c) for r,c in moves if 1 <= r <= 3 and 0 <= c <= 3]

def board_to_string(board):
    return ','.join([''.join(row) for row in board])

def get_position_name(pos):
    row, col = pos
    return f"{chr(65+col)}{row+1}"

def solve_knight_swap():
    # Initial board setup (0-based indexing)
    initial_board = [
        ['.', '.', 'w', '.'],  # row 3 (index 0)
        ['.', '.', '.', 'w'],  # row 2 (index 1)
        ['.', 'B', '.', 'B']   # row 1 (index 2)
    ]
    
    # Target positions (using 0-based indexing)
    white_targets = {(0,2): (2,1), (1,3): (2,3)}  # current_pos: target_pos
    black_targets = {(2,1): (0,2), (2,3): (1,3)}  # current_pos: target_pos
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
        visited.add(board_state)
        
        # Check if solved
        solved = True
        for pos, target in white_targets.items():
            if board[pos[0]][pos[1]] == 'w' or board[target[0]][target[1]] != 'w':
                solved = False
                break
        if solved:
            return moves
        
        # Generate moves
        piece = 'w' if is_white_turn else 'B'
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    curr_pos = (r, c)
                    for new_r, new_c in get_knight_moves((r+1, c)):  # Convert to 1-based for move generation
                        if 0 <= new_r < 3 and 0 <= new_c < 4 and board[new_r][new_c] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[new_r][new_c] = piece
                            new_board[r][c] = '.'
                            move = f"{piece},{get_position_name((r,c))},{get_position_name((new_r,new_c))}"
                            queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")