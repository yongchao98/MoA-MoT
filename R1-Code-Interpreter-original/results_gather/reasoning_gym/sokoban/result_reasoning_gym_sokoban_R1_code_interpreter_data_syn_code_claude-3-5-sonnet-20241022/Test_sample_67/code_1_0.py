def verify_move(sequence):
    # Initial state
    board = [
        list('+++++++++'),
        list('+---------+'),
        list('+-X-X-@--+'),
        list('+--@--$--+'),
        list('+-@*@X---+'),
        list('+X-@$----+'),
        list('+-----$--+'),
        list('+X-------+'),
        list('+++++++++')
    ]
    
    def get_player():
        for i in range(len(board)):
            for j in range(len(board[i])):
                if board[i][j] in ['*', '%']:
                    return (i, j)
        return None
    
    def move(direction):
        player = get_player()
        if not player:
            return False
            
        px, py = player
        dx, dy = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}[direction]
        
        nx, ny = px + dx, py + dy
        
        # Check wall
        if board[nx][ny] == '+':
            return False
            
        # Moving to empty space or goal
        if board[nx][ny] in ['-', 'X']:
            was_on_goal = board[px][py] == '%'
            going_to_goal = board[nx][ny] == 'X'
            
            board[px][py] = 'X' if was_on_goal else '-'
            board[nx][ny] = '%' if going_to_goal else '*'
            return True
            
        # Pushing box
        if board[nx][ny] in ['@', '$']:
            box_x, box_y = nx + dx, ny + dy
            
            if board[box_x][box_y] in ['+', '@', '$']:
                return False
                
            was_on_goal = board[px][py] == '%'
            box_was_on_goal = board[nx][ny] == '$'
            box_going_to_goal = board[box_x][box_y] == 'X'
            
            board[px][py] = 'X' if was_on_goal else '-'
            board[nx][ny] = '%' if box_was_on_goal else '*'
            board[box_x][box_y] = '$' if box_going_to_goal else '@'
            return True
            
        return False
    
    def is_solved():
        boxes = 0
        goals = 0
        boxes_on_goals = 0
        for row in board:
            for cell in row:
                if cell in ['X', '$', '%']:
                    goals += 1
                if cell in ['@', '$']:
                    boxes += 1
                if cell == '$':
                    boxes_on_goals += 1
        return boxes == goals == boxes_on_goals
    
    # Execute moves
    for m in sequence:
        if not move(m):
            return False
    
    return is_solved()

# Test a new solution
solution = "LLUURRDDLLUURRDDLLUURR"
print(solution if verify_move(solution) else "Invalid solution")