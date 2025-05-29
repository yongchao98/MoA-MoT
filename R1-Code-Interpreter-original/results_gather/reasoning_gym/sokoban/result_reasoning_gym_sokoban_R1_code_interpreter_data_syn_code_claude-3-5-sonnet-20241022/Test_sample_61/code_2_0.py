def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_valid_move(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def try_move(board, path, depth, max_depth):
    if depth >= max_depth:
        return None
        
    player = get_player_pos(board)
    if not player:
        return None
        
    px, py = player
    
    # Check if won
    boxes = sum(row.count('@') + row.count('$') for row in board)
    goals_filled = sum(row.count('$') for row in board)
    if boxes == goals_filled and boxes > 0:
        return path
    
    # Try each direction
    directions = [
        ('U', -1, 0),
        ('D', 1, 0),
        ('L', 0, -1),
        ('R', 0, 1)
    ]
    
    for move, dx, dy in directions:
        new_x, new_y = px + dx, py + dy
        
        if not is_valid_move(board, new_x, new_y):
            continue
            
        # Create new board state
        new_board = [list(row) for row in board]
        
        # Moving to empty space or goal
        if new_board[new_x][new_y] in ['-', 'X']:
            new_board[new_x][new_y] = '*' if new_board[new_x][new_y] == '-' else '%'
            new_board[px][py] = '-' if new_board[px][py] == '*' else 'X'
            
            result = try_move(new_board, path + move, depth + 1, max_depth)
            if result:
                return result
                
        # Pushing a box
        elif new_board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if is_valid_move(board, box_x, box_y) and board[box_x][box_y] in ['-', 'X']:
                # Move box
                new_board[box_x][box_y] = '@' if board[box_x][box_y] == '-' else '$'
                # Move player
                new_board[new_x][new_y] = '*' if board[new_x][new_y] == '@' else '%'
                new_board[px][py] = '-' if board[px][py] == '*' else 'X'
                
                result = try_move(new_board, path + move, depth + 1, max_depth)
                if result:
                    return result
    
    return None

# Initial board
board = [
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

# Try with increasing depth limits
for max_depth in [10, 15, 20, 25]:
    solution = try_move(board, "", 0, max_depth)
    if solution:
        print(solution)
        break
    print(f"No solution found at depth {max_depth}")