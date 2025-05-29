from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(board, player_pos):
    rows, cols = len(board), len(board[0])
    moves = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    for dx, dy, move in moves:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        new_board = copy.deepcopy(board)
        new_player_pos = (new_x, new_y)
        
        # If there's a box
        if board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if not is_valid(box_x, box_y, rows, cols) or board[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move the box
            if board[new_x][new_y] == '@':
                new_board[new_x][new_y] = '-'
            else:  # '$'
                new_board[new_x][new_y] = 'X'
                
            if board[box_x][box_y] == 'X':
                new_board[box_x][box_y] = '$'
            else:
                new_board[box_x][box_y] = '@'
        
        # Update player position
        if board[player_pos[0]][player_pos[1]] == '*':
            new_board[player_pos[0]][player_pos[1]] = 'X'
        else:
            new_board[player_pos[0]][player_pos[1]] = '-'
            
        if board[new_x][new_y] == 'X':
            new_board[new_x][new_y] = '%'
        else:
            new_board[new_x][new_y] = '*'
            
        next_states.append((new_board, new_player_pos, move))
    
    return next_states

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's a box not on goal or an empty goal
                return False
    return True

def solve_sokoban(initial_board):
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
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_pos, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        for new_board, new_pos, move in get_next_states(current_board, current_pos):
            board_str = board_to_string(new_board)
            if board_str not in visited:
                visited.add(board_str)
                queue.append((new_board, new_pos, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', 'X', '-', '$', 'X', '+'],
    ['+', '-', '-', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '@', '@', '-', '@', '-', '-', '+'],
    ['+', '-', '@', '-', '-', '@', 'X', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '@', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '*', '+', '-', '$', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")