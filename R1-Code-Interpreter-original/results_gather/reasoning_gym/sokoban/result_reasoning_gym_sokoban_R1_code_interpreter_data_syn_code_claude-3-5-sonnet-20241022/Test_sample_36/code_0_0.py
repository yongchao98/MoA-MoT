from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(board, rows, cols):
    # Find player position
    player_x, player_y = None, None
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['@', '*', '%']:
                player_x, player_y = i, j
                break
        if player_x is not None:
            break

    moves = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    for dx, dy, move in moves:
        new_x, new_y = player_x + dx, player_y + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        new_board = [row[:] for row in board]
        
        # Handle player on goal
        if board[player_x][player_y] == '%':
            new_board[player_x][player_y] = 'X'
        elif board[player_x][player_y] == '*':
            new_board[player_x][player_y] = '-'
        else:
            new_board[player_x][player_y] = '-'
            
        if board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if not is_valid(box_x, box_y, rows, cols) or board[box_x][box_y] in ['+', '@', '$']:
                continue
                
            if board[box_x][box_y] == 'X':
                new_board[box_x][box_y] = '$'
            else:
                new_board[box_x][box_y] = '@'
                
            if board[new_x][new_y] == '$':
                new_board[new_x][new_y] = '%'
            else:
                new_board[new_x][new_y] = '@'
        else:
            if board[new_x][new_y] == 'X':
                new_board[new_x][new_y] = '%'
            else:
                new_board[new_x][new_y] = '@'
                
        next_states.append((new_board, move))
    
    return next_states

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def solve_sokoban(initial_board):
    rows = len(initial_board)
    cols = len(initial_board[0])
    
    queue = deque([(initial_board, "")])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        # Check if solved
        goals_covered = True
        for i in range(rows):
            for j in range(cols):
                if current_board[i][j] == 'X':
                    goals_covered = False
                    break
            if not goals_covered:
                break
                
        if goals_covered:
            return path
            
        for next_board, move in get_next_states(current_board, rows, cols):
            board_str = board_to_string(next_board)
            if board_str not in visited:
                visited.add(board_str)
                queue.append((next_board, path + move))
    
    return None

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '-', '$', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '-', '@', '+', '-', '+'],
    ['+', '$', '$', 'X', '-', '@', '+'],
    ['+', '+', '+', '+', '$', '*', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)