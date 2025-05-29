from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_goal_state(board):
    boxes = 0
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['@', '$']:
                boxes += 1
            if cell in ['X', '$', '%']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return boxes == goals == boxes_on_goals

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    visited = set()
    queue = deque([(initial_board, "", get_player_pos(initial_board))])
    
    while queue:
        board, path, (player_x, player_y) = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
            
        visited.add(board_state)
        
        if is_goal_state(board):
            return path
            
        for dx, dy, move in moves:
            new_x, new_y = player_x + dx, player_y + dy
            
            if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
                continue
                
            new_board = [list(row) for row in board]
            
            if board[new_x][new_y] in ['@', '$']:
                box_x, box_y = new_x + dx, new_y + dy
                
                if not is_valid(box_x, box_y, rows, cols) or \
                   board[box_x][box_y] in ['+', '@', '$']:
                    continue
                    
                # Move box
                new_board[box_x][box_y] = '$' if board[box_x][box_y] in ['X', '%'] else '@'
                new_board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
                new_board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
                
                queue.append((new_board, path + move, (new_x, new_y)))
            else:
                # Move player
                new_board[new_x][new_y] = '%' if board[new_x][new_y] in ['X', '%'] else '*'
                new_board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
                
                queue.append((new_board, path + move, (new_x, new_y)))
    
    return None

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '-', '+', '+', '$', '+', '+', '+'],
    ['+', 'X', '-', '@', '-', '-', 'X', '$', '+'],
    ['+', '@', '-', '+', '$', 'X', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+', 'X', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '$', '-', '-', '-', '+'],
    ['+', '-', '@', '-', '@', '-', '-', '@', '+'],
    ['+', '-', '-', '-', '*', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")