from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at black's initial positions and vice versa
    return (board[0][2] == 'B' and board[1][1] == 'B' and 
            board[2][1] == 'w' and board[0][3] == 'w')

def find_knights(board, color):
    knights = []
    for i in range(4):
        for j in range(3):
            if board[j][i] == color:
                knights.append((i, j))
    return knights

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', '.', 'B', 'w'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(current_board):
            return moves
        
        color = 'w' if is_white_turn else 'B'
        knights = find_knights(current_board, color)
        
        for kx, ky in knights:
            for new_x, new_y in get_knight_moves((kx, ky)):
                if current_board[new_y][new_x] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[ky][kx] = '.'
                    new_board[new_y][new_x] = color
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{color},{chr(65+kx)}{ky+1},{chr(65+new_x)}{new_y+1}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")