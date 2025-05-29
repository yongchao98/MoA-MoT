def verify_solution(board, moves):
    height, width = len(board), len(board[0])
    # Deep copy the board
    current = [row[:] for row in board]
    
    # Find player position
    player_x, player_y = None, None
    for i in range(height):
        for j in range(width):
            if current[i][j] in ['@', '*']:
                player_x, player_y = i, j
                break
        if player_x is not None:
            break
    
    def move(dx, dy):
        nonlocal player_x, player_y
        new_x, new_y = player_x + dx, player_y + dy
        
        # Check if move is valid
        if current[new_x][new_y] == '+':
            return False
            
        # If there's a box
        if current[new_x][new_y] in ['$', '@']:
            next_x, next_y = new_x + dx, new_y + dy
            if current[next_x][next_y] in ['+', '$', '@']:
                return False
            # Move box
            if current[new_x][new_y] == '$':
                current[next_x][next_y] = '$'
            current[new_x][new_y] = '@'
        else:
            current[new_x][new_y] = '@'
            
        # Update old position
        current[player_x][player_y] = '-'
        player_x, player_y = new_x, new_y
        return True
    
    # Try the solution
    for m in moves:
        if m == 'U':
            if not move(-1, 0): return False
        elif m == 'D':
            if not move(1, 0): return False
        elif m == 'L':
            if not move(0, -1): return False
        elif m == 'R':
            if not move(0, 1): return False
    
    # Check if all boxes are on goals
    goals = 0
    boxes = 0
    for i in range(height):
        for j in range(width):
            if current[i][j] == '$':
                boxes += 1
            if board[i][j] == 'X':
                goals += 1
                if current[i][j] != '$':
                    return False
    
    return boxes == goals

# Test board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '$', '$', '-', '@', '-', '+'],
    ['+', '-', '-', '-', 'X', 'X', 'X', '-', '-', '+'],
    ['+', '+', 'X', '-', '@', '@', '@', '@', '*', '+'],
    ['+', '+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

# Test solution
solution = "LLLUURDDLLUULLDDRR"
print(verify_solution(board, solution))