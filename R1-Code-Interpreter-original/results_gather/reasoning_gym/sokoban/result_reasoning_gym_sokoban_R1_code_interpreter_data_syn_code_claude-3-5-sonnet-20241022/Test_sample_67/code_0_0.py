from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '@']:
                return (i, j)
    return None

def is_goal_complete(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['X', '$', '%']:
                goals += 1
            if cell in ['$']:
                boxes_on_goals += 1
    return goals == boxes_on_goals

def get_next_states(board, visited):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    states = []
    player = get_player_pos(board)
    if not player:
        return states
    
    px, py = player
    
    for dx, dy, move in directions:
        new_px, new_py = px + dx, py + dy
        if not is_valid(new_px, new_py, len(board), len(board[0])):
            continue
            
        if board[new_px][new_py] in ['@', '$']:  # pushing a box
            box_new_x, box_new_y = new_px + dx, new_py + dy
            if (is_valid(box_new_x, box_new_y, len(board), len(board[0])) and
                board[box_new_x][box_new_y] not in ['+', '@', '$']):
                
                new_board = [list(row) for row in board]
                # Move box
                if board[new_px][new_py] == '@':
                    new_board[box_new_x][box_new_y] = '$' if board[box_new_x][box_new_y] in ['X'] else '@'
                else:  # box was on goal
                    new_board[box_new_x][box_new_y] = '$' if board[box_new_x][box_new_y] in ['X'] else '@'
                
                # Move player
                new_board[new_px][new_py] = '%' if board[new_px][new_py] in ['$'] else '*'
                new_board[px][py] = 'X' if board[px][py] in ['%'] else '-'
                
                board_tuple = tuple(map(tuple, new_board))
                if board_tuple not in visited:
                    states.append((new_board, move))
                    
        elif board[new_px][new_py] not in ['+']:  # moving to empty space or goal
            new_board = [list(row) for row in board]
            new_board[new_px][new_py] = '%' if board[new_px][new_py] in ['X'] else '*'
            new_board[px][py] = 'X' if board[px][py] in ['%'] else '-'
            
            board_tuple = tuple(map(tuple, new_board))
            if board_tuple not in visited:
                states.append((new_board, move))
    
    return states

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    visited = {tuple(map(tuple, initial_board))}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_goal_complete(current_board):
            return moves
            
        for next_board, move in get_next_states(current_board, visited):
            visited.add(tuple(map(tuple, next_board)))
            queue.append((next_board, moves + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', 'X', '-', '@', '-', '+'],
    ['+', '-', '-', '@', '-', '-', '$', '-', '+'],
    ['+', '-', '@', '*', '@', 'X', '-', '-', '+'],
    ['+', 'X', '-', '@', '$', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '$', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)