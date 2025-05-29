from collections import deque

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:  # Player or player on goal
                return (i, j)
    return None

def get_boxes(board):
    boxes = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:  # Box or box on goal
                boxes.append((i, j))
    return boxes

def get_goals(board):
    goals = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['X', '$', '%']:  # Goal, box on goal, or player on goal
                goals.append((i, j))
    return goals

def is_valid_pos(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def make_move(board, player_pos, dx, dy):
    px, py = player_pos
    new_x, new_y = px + dx, py + dy
    
    if not is_valid_pos(board, new_x, new_y):
        return None
        
    new_board = [row[:] for row in board]
    
    # Check if moving to empty space or goal
    if board[new_x][new_y] in ['-', 'X']:
        # Update player position
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[new_x][new_y] = '%' if board[new_x][new_y] == 'X' else '*'
        return new_board
        
    # Check if pushing a box
    if board[new_x][new_y] in ['@', '$']:
        push_x, push_y = new_x + dx, new_y + dy
        if is_valid_pos(board, push_x, push_y) and board[push_x][push_y] in ['-', 'X']:
            # Update player and box positions
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            new_board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
            new_board[push_x][push_y] = '$' if board[push_x][push_y] == 'X' else '@'
            return new_board
            
    return None

def is_solved(board):
    goals = set(map(tuple, get_goals(board)))
    boxes = set(map(tuple, get_boxes(board)))
    return all(box in goals for box in boxes)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(board):
    directions = [
        (0, 1, 'R'),   # Right
        (0, -1, 'L'),  # Left
        (-1, 0, 'U'),  # Up
        (1, 0, 'D')    # Down
    ]
    
    initial_state = [row[:] for row in board]
    queue = deque([(initial_state, "")])
    visited = {board_to_string(initial_state)}
    max_moves = 200  # Limit search depth
    
    while queue:
        current_board, moves = queue.popleft()
        
        if len(moves) >= max_moves:
            continue
            
        if is_solved(current_board):
            return moves
            
        player_pos = get_player_pos(current_board)
        if not player_pos:
            continue
            
        for dx, dy, direction in directions:
            new_board = make_move(current_board, player_pos, dx, dy)
            if new_board:
                board_str = board_to_string(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    queue.append((new_board, moves + direction))
    
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