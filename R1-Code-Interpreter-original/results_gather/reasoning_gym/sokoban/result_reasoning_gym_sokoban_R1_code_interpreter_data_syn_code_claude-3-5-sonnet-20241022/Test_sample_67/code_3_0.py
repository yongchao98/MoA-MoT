def verify_solution(moves):
    # Initial board state
    board = [
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
    
    def print_board():
        for row in board:
            print(''.join(row))
    
    def get_player_pos():
        for i in range(len(board)):
            for j in range(len(board[0])):
                if board[i][j] in ['*', '%']:
                    return (i, j)
        return None
    
    def move(direction):
        px, py = get_player_pos()
        dx, dy = {
            'U': (-1, 0),
            'D': (1, 0),
            'L': (0, -1),
            'R': (0, 1)
        }[direction]
        
        nx, ny = px + dx, py + dy
        
        # Check if move is valid
        if board[nx][ny] == '+':
            return False
            
        if board[nx][ny] in ['@', '$']:  # pushing a box
            box_x, box_y = nx + dx, ny + dy
            if board[box_x][box_y] in ['+', '@', '$']:
                return False
                
            # Move box
            is_goal = board[box_x][box_y] == 'X'
            board[box_x][box_y] = '$' if is_goal else '@'
            
            # Move player
            is_from_goal = board[px][py] == '%'
            board[px][py] = 'X' if is_from_goal else '-'
            board[nx][ny] = '%' if board[nx][ny] == '$' else '*'
            
        else:  # moving to empty space
            is_from_goal = board[px][py] == '%'
            is_to_goal = board[nx][ny] == 'X'
            board[px][py] = 'X' if is_from_goal else '-'
            board[nx][ny] = '%' if is_to_goal else '*'
            
        return True
    
    def is_solved():
        goals = 0
        boxes_on_goals = 0
        for row in board:
            for cell in row:
                if cell in ['X', '$', '%']:
                    goals += 1
                if cell == '$':
                    boxes_on_goals += 1
        return goals == boxes_on_goals
    
    # Execute moves
    for m in moves:
        if not move(m):
            return False
    
    return is_solved()

# Test solution
solution = "ULLDDRRUULLDDRRUULLDDR"
print(solution if verify_solution(solution) else "Invalid solution")