def is_valid_move(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def can_push_box(board, player_x, player_y, box_x, box_y, dx, dy):
    new_box_x, new_box_y = box_x + dx, box_y + dy
    return is_valid_move(board, new_box_x, new_box_y) and board[new_box_x][new_box_y] in ['-', 'X']

def try_solution(moves):
    board = [
        list("++++++"),
        list("+*@-X+"),
        list("+@@$++"),
        list("+X-$++"),
        list("++-+$+"),
        list("+--+++"),
        list("+--X++"),
        list("+-++++"),
        list("++++++"),
        list("++++++")
    ]
    
    player_x, player_y = 1, 1  # Initial player position
    
    for move in moves:
        dx, dy = 0, 0
        if move == 'U': dx = -1
        elif move == 'D': dx = 1
        elif move == 'L': dy = -1
        elif move == 'R': dy = 1
        
        new_x, new_y = player_x + dx, player_y + dy
        
        if not is_valid_move(board, new_x, new_y):
            return False
            
        if board[new_x][new_y] in ['@', '$']:  # If there's a box
            if not can_push_box(board, player_x, player_y, new_x, new_y, dx, dy):
                return False
                
            # Move the box
            box_new_x, box_new_y = new_x + dx, new_y + dy
            is_goal = board[box_new_x][box_new_y] == 'X'
            board[box_new_x][box_new_y] = '$' if is_goal else '@'
            
            # Update previous box position
            board[new_x][new_y] = '*'
            
            # Update previous player position
            board[player_x][player_y] = '-' if board[player_x][player_y] == '*' else 'X'
            
        else:  # Simple player move
            board[new_x][new_y] = '*'
            board[player_x][player_y] = '-' if board[player_x][player_y] == '*' else 'X'
            
        player_x, player_y = new_x, new_y
        
    # Check if all boxes are on goals
    for row in board:
        for cell in row:
            if cell in ['@', 'X']:  # If there's a box not on goal or an empty goal
                return False
    return True

# Test the specific solution
solution = "DLDRRUULLDDRRUULDDLLUURRDDLLUURR"
if try_solution(solution):
    print(solution)
else:
    print("Solution not valid")