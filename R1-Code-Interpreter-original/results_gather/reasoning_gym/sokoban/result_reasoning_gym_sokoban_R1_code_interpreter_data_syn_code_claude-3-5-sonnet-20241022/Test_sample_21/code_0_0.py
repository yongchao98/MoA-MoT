from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_state_string(board):
    return '\n'.join([''.join(row) for row in board])

def is_goal(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's a box not on goal or goal without box
                return False
    return True

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    queue = deque([(initial_board, "", get_player_pos(initial_board))])
    visited = set()
    
    while queue:
        board, path, (player_x, player_y) = queue.popleft()
        state = get_state_string(board)
        
        if state in visited:
            continue
        
        visited.add(state)
        
        if is_goal(board):
            return path
        
        for dx, dy, move in moves:
            new_x, new_y = player_x + dx, player_y + dy
            
            if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
                continue
            
            new_board = [list(row) for row in board]
            
            if board[new_x][new_y] in ['@', '$']:  # If pushing a box
                box_x, box_y = new_x + dx, new_y + dy
                
                if not is_valid(box_x, box_y, rows, cols) or \
                   board[box_x][box_y] in ['+', '@', '$']:
                    continue
                
                # Move box
                new_board[box_x][box_y] = '$' if board[box_x][box_y] == 'X' else '@'
                new_board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
                new_board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
                
                queue.append((new_board, path + move, (new_x, new_y)))
            else:  # Just moving player
                new_board[new_x][new_y] = '%' if board[new_x][new_y] == 'X' else '*'
                new_board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
                
                queue.append((new_board, path + move, (new_x, new_y)))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '*', '@', 'X', '+', '+', '-', '-', '+', '+'],
    ['+', '$', '-', '+', '+', '-', '@', '-', '+', '+'],
    ['+', '$', '-', '$', 'X', '-', '-', '-', '+', '+'],
    ['+', '$', '-', '@', '-', '-', 'X', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)