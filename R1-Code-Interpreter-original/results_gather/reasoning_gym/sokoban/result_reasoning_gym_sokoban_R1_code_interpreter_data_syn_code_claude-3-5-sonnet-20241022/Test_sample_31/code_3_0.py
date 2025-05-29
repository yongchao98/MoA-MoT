def verify_move(board, move_sequence):
    # Convert board to list of lists for mutability
    board = [list(row) for row in board]
    
    # Find player position
    player = None
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == '*':
                player = [i, j]
                break
        if player:
            break
    
    # Define move directions
    moves = {
        'U': [-1, 0],
        'D': [1, 0],
        'L': [0, -1],
        'R': [0, 1]
    }
    
    # Try each move
    for move in move_sequence:
        if move not in moves:
            return False
            
        dx, dy = moves[move]
        new_x, new_y = player[0] + dx, player[1] + dy
        
        # Check wall collision
        if board[new_x][new_y] == '+':
            return False
            
        # Moving to empty space
        if board[new_x][new_y] == '-':
            board[player[0]][player[1]] = '-'
            board[new_x][new_y] = '*'
            player = [new_x, new_y]
            continue
            
        # Moving to goal
        if board[new_x][new_y] == 'X':
            board[player[0]][player[1]] = '-'
            board[new_x][new_y] = '*'
            player = [new_x, new_y]
            continue
            
        # Pushing box
        if board[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            
            # Check if box can be pushed
            if board[push_x][push_y] in ['+', '@', '$']:
                return False
                
            # Push box
            if board[push_x][push_y] == 'X':
                board[push_x][push_y] = '$'
            else:
                board[push_x][push_y] = '@'
                
            board[new_x][new_y] = '*'
            board[player[0]][player[1]] = '-'
            player = [new_x, new_y]
            
    return True

# Test board
board = [
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

# Let's try a shorter solution
solution = "ULUURRDDLLUURR"
print(solution if verify_move(board, solution) else "Invalid solution")