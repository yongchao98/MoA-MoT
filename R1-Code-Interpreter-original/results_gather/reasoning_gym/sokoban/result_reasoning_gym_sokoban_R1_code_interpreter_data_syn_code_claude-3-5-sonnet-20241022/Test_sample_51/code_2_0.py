from collections import deque

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '*', '%']:
                return (i, j)
    return None

def is_valid(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def is_complete(board):
    goals = sum(row.count('X') + row.count('$') + row.count('%') for row in board)
    boxes_on_goals = sum(row.count('$') + row.count('%') for row in board)
    return goals == boxes_on_goals

def make_move(board, player_pos, dx, dy):
    px, py = player_pos
    new_x, new_y = px + dx, py + dy
    
    if not is_valid(board, new_x, new_y):
        return None
        
    new_board = [list(row) for row in board]
    
    # Moving to empty space or goal
    if board[new_x][new_y] in ['-', 'X']:
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[new_x][new_y] = '%' if board[new_x][new_y] == 'X' else '*'
        return new_board
        
    # Pushing a box
    if board[new_x][new_y] in ['@', '$']:
        push_x, push_y = new_x + dx, new_y + dy
        if is_valid(board, push_x, push_y) and board[push_x][push_y] in ['-', 'X']:
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            new_board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
            new_board[push_x][push_y] = '$' if board[push_x][push_y] == 'X' else '@'
            return new_board
            
    return None

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(board):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    initial = [row[:] for row in board]
    queue = deque([(initial, "")])
    visited = {board_to_string(initial)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_complete(current_board):
            return path
            
        player_pos = get_player_pos(current_board)
        if not player_pos:
            continue
            
        for dx, dy, move in directions:
            new_board = make_move(current_board, player_pos, dx, dy)
            if new_board:
                board_str = board_to_string(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    queue.append((new_board, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '+', 'X', '+', '$', 'X', '+'],
    ['+', '-', 'X', 'X', '@', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '-', '-', '+'],
    ['+', 'X', '-', '@', '@', '@', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '@', '@', '-', '-', '-', '+', '+'],
    ['+', '-', '%', '-', '-', '-', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")