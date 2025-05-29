from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row-2, col+1), (row-2, col-1),
        (row+2, col+1), (row+2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def pos_to_notation(pos):
    row, col = pos
    return chr(65 + col) + str(3 - row)

def state_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at black's initial positions and vice versa
    return (board[0][2] == 'w' and board[2][2] == 'w' and 
            board[1][1] == 'B' and board[2][3] == 'B')

def solve_knight_swap():
    # Initial board state (0,0 is top-left)
    initial_board = [
        ['.', '.', 'B', '.'],  # row 3
        ['.', 'w', '.', '.'],  # row 2
        ['.', '.', 'B', 'w']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {state_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find current positions of the piece
        positions = [(r, c) for r in range(3) for c in range(4) if board[r][c] == piece]
        
        # Try all possible moves for each piece
        for pos in positions:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece
                    
                    board_string = state_to_string(new_board)
                    if board_string not in visited:
                        visited.add(board_string)
                        move = f"{piece},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")