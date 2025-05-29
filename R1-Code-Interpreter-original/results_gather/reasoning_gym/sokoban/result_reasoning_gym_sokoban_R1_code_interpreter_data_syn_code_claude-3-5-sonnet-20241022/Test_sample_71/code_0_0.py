from collections import deque
import copy

def is_valid_move(board, player_pos, direction):
    rows, cols = len(board), len(board[0])
    r, c = player_pos
    dr, dc = direction
    new_r, new_c = r + dr, c + dc
    
    # Check if new position is within bounds and not a wall
    if not (0 <= new_r < rows and 0 <= new_c < cols) or board[new_r][new_c] == '+':
        return False
        
    # If moving to empty space or goal
    if board[new_r][new_c] in ['-', 'X']:
        return True
        
    # If pushing a box, check if space behind box is available
    if board[new_r][new_c] in ['@', '$']:
        next_r, next_c = new_r + dr, new_c + dc
        if not (0 <= next_r < rows and 0 <= next_c < cols):
            return False
        return board[next_r][next_c] in ['-', 'X']
        
    return False

def make_move(board, player_pos, direction):
    new_board = [list(row) for row in board]
    r, c = player_pos
    dr, dc = direction
    new_r, new_c = r + dr, c + dc
    
    # Remove player from current position
    if board[r][c] == '*':
        new_board[r][c] = 'X'
    else:
        new_board[r][c] = '-'
    
    # If pushing a box
    if board[new_r][new_c] in ['@', '$']:
        next_r, next_c = new_r + dr, new_c + dc
        # Move box
        if board[next_r][next_c] == 'X':
            new_board[next_r][next_c] = '$'
        else:
            new_board[next_r][next_c] = '@'
        
    # Move player
    if board[new_r][new_c] in ['@', '-']:
        new_board[new_r][new_c] = '*'
    else:  # moving to a goal
        new_board[new_r][new_c] = '%'
        
    return new_board, (new_r, new_c)

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's a box not on goal or empty goal
                return False
    return True

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    # Find initial player position
    player_pos = None
    for i in range(len(initial_board)):
        for j in range(len(initial_board[0])):
            if initial_board[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_board, player_pos, "")])
    visited = set()
    
    while queue:
        current_board, pos, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        board_state = board_to_string(current_board)
        if board_state in visited:
            continue
        visited.add(board_state)
        
        for move, dir_vec in directions.items():
            if is_valid_move(current_board, pos, dir_vec):
                new_board, new_pos = make_move(current_board, pos, dir_vec)
                queue.append((new_board, new_pos, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '+'],
    ['+', 'X', '@', '-', 'X', '-', '+'],
    ['+', '-', '@', '*', '@', '-', '+'],
    ['+', '-', '-', 'X', '@', 'X', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)