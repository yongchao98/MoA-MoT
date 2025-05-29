def is_valid_move(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def try_solution(moves):
    # Initialize the board
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
    
    # For debugging
    print("Initial state:")
    for row in board:
        print(''.join(row))
    print()
    
    for i, move in enumerate(moves):
        dx, dy = 0, 0
        if move == 'U': dx = -1
        elif move == 'D': dx = 1
        elif move == 'L': dy = -1
        elif move == 'R': dy = 1
        
        new_x, new_y = player_x + dx, player_y + dy
        
        # Check if move is valid
        if not is_valid_move(board, new_x, new_y):
            print(f"Invalid move at step {i+1}: {move}")
            return False
            
        # Moving into a box
        if board[new_x][new_y] in ['@', '$']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            
            if not is_valid_move(board, box_new_x, box_new_y):
                print(f"Cannot push box at step {i+1}: {move}")
                return False
            
            # Move box
            was_on_goal = board[new_x][new_y] == '$'
            moving_to_goal = board[box_new_x][box_new_y] == 'X'
            board[box_new_x][box_new_y] = '$' if moving_to_goal else '@'
            
            # Move player
            board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
            board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
            
        else:
            # Simple player move
            board[new_x][new_y] = '%' if board[new_x][new_y] == 'X' else '*'
            board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
        
        player_x, player_y = new_x, new_y
        
        # Debug: print state after each move
        print(f"After move {i+1}: {move}")
        for row in board:
            print(''.join(row))
        print()
    
    # Check win condition
    for row in board:
        for cell in row:
            if cell in ['@', 'X']:  # If there's a box not on goal or an empty goal
                return False
    return True

# Test a new solution
solution = "RDLDDRRULLDDRRULLDDRRUULL"
if try_solution(solution):
    print("Solution valid!")
else:
    print("Solution not valid")