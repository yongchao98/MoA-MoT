def print_board(board):
    for row in board:
        print(row)
    print()

def make_move(board, direction):
    board = [list(row) for row in board]
    
    # Find player
    player_pos = None
    for y in range(len(board)):
        for x in range(len(board[y])):
            if board[y][x] in ['*', '%']:
                player_pos = (x, y)
                break
        if player_pos:
            break
    
    if not player_pos:
        return None
    
    px, py = player_pos
    
    # Define movement directions
    if direction == 'L':
        dx, dy = -1, 0
    elif direction == 'R':
        dx, dy = 1, 0
    elif direction == 'U':
        dx, dy = 0, -1
    elif direction == 'D':
        dx, dy = 0, 1
    else:
        return None
    
    # Calculate new positions
    new_px = px + dx
    new_py = py + dy
    
    # Check if move is valid
    if board[new_py][new_px] == '+':
        return None
    
    # Moving to empty space or goal
    if board[new_py][new_px] in ['-', 'X']:
        # Update player's old position
        board[py][px] = '-' if board[py][px] == '*' else 'X'
        # Update player's new position
        board[new_py][new_px] = '*' if board[new_py][new_px] == '-' else '%'
    
    # Pushing a box
    elif board[new_py][new_px] in ['@', '$']:
        box_x = new_px + dx
        box_y = new_py + dy
        
        # Check if box can be pushed
        if box_x < 0 or box_y < 0 or box_x >= len(board[0]) or box_y >= len(board):
            return None
        if board[box_y][box_x] in ['+', '@', '$']:
            return None
            
        # Move the box
        is_box_on_goal = board[new_py][new_px] == '$'
        is_new_pos_goal = board[box_y][box_x] == 'X'
        
        # Update player's old position
        board[py][px] = '-' if board[py][px] == '*' else 'X'
        
        # Update player's new position
        board[new_py][new_px] = '*' if not is_box_on_goal else '%'
        
        # Update box's new position
        board[box_y][box_x] = '@' if not is_new_pos_goal else '$'
    
    return [''.join(row) for row in board]

def verify_solution(solution):
    board = [
        "+++++++",
        "++X@*+",
        "++XX@-+",
        "++@-+",
        "++$+-+",
        "+++++++"
    ]
    
    print("Initial state:")
    print_board(board)
    
    for move in solution:
        print(f"Making move: {move}")
        new_board = make_move(board, move)
        if not new_board:
            print(f"Invalid move: {move}")
            return False
        board = new_board
        print_board(board)
    
    # Check if all boxes are on goals
    boxes_on_goals = 0
    total_goals = 0
    for row in board:
        for cell in row:
            if cell in ['X', '%']:
                total_goals += 1
            if cell in ['$', '%']:
                boxes_on_goals += 1
    
    return boxes_on_goals == total_goals and boxes_on_goals == 3

# Test different solutions
solutions_to_try = [
    "LULDRRL",
    "LUURDLL",
    "LLURRDL"
]

for solution in solutions_to_try:
    print(f"\nTrying solution: {solution}")
    if verify_solution(solution):
        print(f"Solution found: {solution}")
        print("<<<" + solution + ">>>")
        break