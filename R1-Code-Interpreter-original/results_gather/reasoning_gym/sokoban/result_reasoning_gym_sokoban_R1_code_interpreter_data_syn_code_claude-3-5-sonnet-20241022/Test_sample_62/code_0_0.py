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

def is_complete(board):
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

def get_next_states(board, visited):
    player = get_player_pos(board)
    if not player:
        return []
    
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    for dx, dy, move in directions:
        new_x, new_y = player[0] + dx, player[1] + dy
        
        if not is_valid(new_x, new_y, len(board), len(board[0])):
            continue
            
        if board[new_x][new_y] == '+':
            continue
            
        new_board = [list(row) for row in board]
        
        # If moving to empty space or goal
        if board[new_x][new_y] in ['-', 'X']:
            new_board[player[0]][player[1]] = '-' if board[player[0]][player[1]] == '*' else 'X'
            new_board[new_x][new_y] = '*' if board[new_x][new_y] == '-' else '%'
            board_tuple = tuple(map(tuple, new_board))
            if board_tuple not in visited:
                next_states.append((new_board, move))
                
        # If pushing a box
        elif board[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            
            if not is_valid(push_x, push_y, len(board), len(board[0])):
                continue
                
            if board[push_x][push_y] in ['-', 'X']:
                new_board[player[0]][player[1]] = '-' if board[player[0]][player[1]] == '*' else 'X'
                new_board[new_x][new_y] = '*' if board[new_x][new_y] == '@' else '%'
                new_board[push_x][push_y] = '@' if board[push_x][push_y] == '-' else '$'
                board_tuple = tuple(map(tuple, new_board))
                if board_tuple not in visited:
                    next_states.append((new_board, move))
                    
    return next_states

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    visited = {tuple(map(tuple, initial_board))}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_complete(current_board):
            return moves
            
        for next_board, move in get_next_states(current_board, visited):
            visited.add(tuple(map(tuple, next_board)))
            queue.append((next_board, moves + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '*', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '@', '@', '@', '@', '-', '-', '@', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '-', '$', '-', '+'],
    ['+', 'X', '-', '-', '-', '$', 'X', '-', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)