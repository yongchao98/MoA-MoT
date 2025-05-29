def print_board_state(board, moves):
    board = [list(row) for row in board]
    player = None
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == '*':
                player = [i, j]
                break
        if player:
            break
    
    directions = {
        'U': [-1, 0],
        'D': [1, 0],
        'L': [0, -1],
        'R': [0, 1]
    }
    
    print("Initial state:")
    for row in board:
        print(''.join(row))
    
    for move in moves:
        dx, dy = directions[move]
        new_x, new_y = player[0] + dx, player[1] + dy
        print(f"\nAfter move {move}:")
        
        # Update board
        if board[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            # Move box
            board[push_x][push_y] = '@'
            board[new_x][new_y] = '*'
            board[player[0]][player[1]] = '-'
        else:
            board[new_x][new_y] = '*'
            board[player[0]][player[1]] = '-'
        
        player = [new_x, new_y]
        
        # Print current state
        for row in board:
            print(''.join(row))

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

# Let's try this solution
test_moves = "ULLUURRDLL"
print_board_state(board, test_moves)