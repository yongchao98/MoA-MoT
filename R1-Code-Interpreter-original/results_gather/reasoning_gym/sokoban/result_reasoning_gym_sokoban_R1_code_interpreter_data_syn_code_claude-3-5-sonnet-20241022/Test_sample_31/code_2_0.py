def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def try_move(board, move_sequence):
    board = [list(row) for row in board]
    player = get_player_pos(board)
    if not player:
        return False
    
    for move in move_sequence:
        px, py = player
        if move == 'U': dx, dy = -1, 0
        elif move == 'D': dx, dy = 1, 0
        elif move == 'L': dx, dy = 0, -1
        elif move == 'R': dx, dy = 0, 1
        
        new_x, new_y = px + dx, py + dy
        
        # Check if move is valid
        if board[new_x][new_y] == '+':
            return False
            
        # Moving to empty space or goal
        if board[new_x][new_y] in ['-', 'X']:
            # Update player position
            if board[px][py] == '*':
                board[px][py] = '-'
            else:  # '%'
                board[px][py] = 'X'
                
            if board[new_x][new_y] == 'X':
                board[new_x][new_y] = '%'
            else:
                board[new_x][new_y] = '*'
            
        # Pushing a box
        elif board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            # Check if box can be pushed
            if board[box_x][box_y] in ['+', '@', '$']:
                return False
                
            # Push the box
            if board[box_x][box_y] == 'X':
                board[box_x][box_y] = '$'
            else:
                board[box_x][box_y] = '@'
                
            if board[new_x][new_y] == '@':
                board[new_x][new_y] = '*'
            else:  # '$'
                board[new_x][new_y] = '%'
                
            if board[px][py] == '*':
                board[px][py] = '-'
            else:  # '%'
                board[px][py] = 'X'
        
        player = (new_x, new_y)
    
    return True

# Initialize the board
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

# Known solution to test
solution = "ULUURRDLULLDDRRUULLDDRRUULLDDRR"
result = try_move(board, solution)
print(solution if result else "Invalid solution")