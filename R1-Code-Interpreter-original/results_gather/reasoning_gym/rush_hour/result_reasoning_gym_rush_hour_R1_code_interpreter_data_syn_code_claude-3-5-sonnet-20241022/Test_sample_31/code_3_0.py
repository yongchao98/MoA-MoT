def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i,j)]
                else:
                    cars[board[i][j]].append((i,j))
    return cars

def make_move(board, car, positions, direction):
    new_board = [list(row) for row in board]
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Determine if car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Add new positions
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = car
    
    return [''.join(row) for row in new_board]

def can_move(board, positions, direction):
    is_horizontal = positions[0][0] == positions[1][0]
    if is_horizontal:
        row = positions[0][0]
        if direction > 0:  # moving right
            rightmost = max(p[1] for p in positions)
            return rightmost + 1 < 6 and board[row][rightmost + 1] == '.'
        else:  # moving left
            leftmost = min(p[1] for p in positions)
            return leftmost > 0 and board[row][leftmost - 1] == '.'
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # moving down
            bottommost = max(p[0] for p in positions)
            return bottommost + 1 < 6 and board[bottommost + 1][col] == '.'
        else:  # moving up
            topmost = min(p[0] for p in positions)
            return topmost > 0 and board[topmost - 1][col] == '.'

def solve_puzzle():
    initial_board = [
        "..xBBL",
        "..ICCL",
        "AAIJ..",
        "HDDJEE",
        "HFFJKx",
        "GG..K."
    ]
    
    # Predetermined solution sequence
    moves = [
        ('J', -1),  # Move J up
        ('K', -1),  # Move K up
        ('E', -1),  # Move E left
        ('C', -1),  # Move C left
        ('I', -1),  # Move I up
        ('L', -1),  # Move L up
        ('B', -1),  # Move B left
        ('A', 4)    # Move A right to exit
    ]
    
    current_board = initial_board
    solution = []
    
    for car, direction in moves:
        cars = get_car_positions(current_board)
        if car in cars and can_move(current_board, cars[car], direction):
            current_board = make_move(current_board, car, cars[car], direction)
            move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
            solution.append(move_str)
        else:
            return None
    
    return solution

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")