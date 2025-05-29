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
        
        if not is_valid_move(board, new_x, new_y):
            print(f"Invalid move at step {i+1}: {move}")
            return False
            
        if board[new_x][new_y] in ['@', '$']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            
            if not is_valid_move(board, box_new_x, box_new_y) or board[box_new_x][box_new_y] in ['@', '$']:
                print(f"Cannot push box at step {i+1}: {move}")
                return False
            
            # Move box
            board[box_new_x][box_new_y] = '$' if board[box_new_x][box_new_y] == 'X' else '@'
            board[new_x][new_y] = '*'
            board[player_x][player_y] = '-'
            
        else:
            # Simple player move
            board[new_x][new_y] = '*'
            board[player_x][player_y] = '-'
        
        player_x, player_y = new_x, new_y
        
        print(f"After move {i+1}: {move}")
        for row in board:
            print(''.join(row))
        print()
    
    # Check win condition
    for row in board:
        for cell in row:
            if cell in ['@', 'X']:
                return False
    return True

# Test a new solution
solution = "RDDLLUURRDDLLUURRDD"
if try_solution(solution):
    print("Solution valid!")
else:
    print("Solution not valid")