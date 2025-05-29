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
    return [(r,c) for r,c in moves if 1 <= r <= 3 and 0 <= c <= 3]

def board_to_string(board):
    return ','.join([''.join(row) for row in board])

def get_position_name(pos):
    row, col = pos
    return f"{chr(65+col)}{row}"

def solve_knight_swap():
    # Initial board setup (1-based rows, 0-based columns)
    # Correct initial state:
    # 3 | . | . | w | . |
    # 2 | . | . | . | w |
    # 1 | . | B | . | B |
    initial_board = [
        ['.' for _ in range(4)] for _ in range(4)]  # Using 0-3 indices
    # Place white knights
    initial_board[3-1][2] = 'w'  # C3
    initial_board[2-1][3] = 'w'  # D2
    # Place black knights
    initial_board[1-1][1] = 'B'  # B1
    initial_board[1-1][3] = 'B'  # D1
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
        visited.add(board_state)
        
        # Check if solved - white knights should be at B1 and D1, black at C3 and D2
        if (board[1-1][1] == 'w' and board[1-1][3] == 'w' and 
            board[3-1][2] == 'B' and board[2-1][3] == 'B'):
            return moves
        
        # Generate moves
        piece = 'w' if is_white_turn else 'B'
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    for new_r, new_c in get_knight_moves((r+1,c)):
                        new_r -= 1  # Convert to 0-based index
                        if board[new_r][new_c] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[new_r][new_c] = piece
                            new_board[r][c] = '.'
                            move = f"{piece},{get_position_name((r+1,c))},{get_position_name((new_r+1,new_c))}"
                            queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")