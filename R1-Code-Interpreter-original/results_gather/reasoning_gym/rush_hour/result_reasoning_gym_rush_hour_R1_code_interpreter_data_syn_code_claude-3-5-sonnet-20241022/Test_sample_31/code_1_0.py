def print_board(board):
    for row in board:
        print(row)

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

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def make_move(board, car, positions, direction):
    new_board = [list(row) for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    for i, j in positions:
        new_i = i if horizontal else i + direction
        new_j = j + direction if horizontal else j
        if 0 <= new_i < 6 and 0 <= new_j < 6:
            new_board[new_i][new_j] = car
        else:
            return None
    
    return [''.join(row) for row in new_board]

def can_move(board, car, positions, direction):
    horizontal = is_horizontal(positions)
    
    if horizontal:
        row = positions[0][0]
        if direction > 0:  # moving right
            rightmost = max(p[1] for p in positions)
            return (rightmost + direction) < 6 and all(board[row][rightmost + i] == '.' for i in range(1, direction + 1))
        else:  # moving left
            leftmost = min(p[1] for p in positions)
            return (leftmost + direction) >= 0 and all(board[row][leftmost + i] == '.' for i in range(direction, 0))
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # moving down
            bottommost = max(p[0] for p in positions)
            return (bottommost + direction) < 6 and all(board[bottommost + i][col] == '.' for i in range(1, direction + 1))
        else:  # moving up
            topmost = min(p[0] for p in positions)
            return (topmost + direction) >= 0 and all(board[topmost + i][col] == '.' for i in range(direction, 0))

def solve_puzzle():
    initial_board = [
        "..xBBL",
        "..ICCL",
        "AAIJ..",
        "HDDJEE",
        "HFFJKx",
        "GG..K."
    ]
    
    # Define the solution sequence with exact steps
    solution_sequence = [
        ('J', -1),  # J up 1
        ('K', -1),  # K up 1
        ('E', -1),  # E left 1
        ('C', -1),  # C left 1
        ('I', -1),  # I up 1
        ('L', -1),  # L up 1
        ('B', -1),  # B left 1
        ('A', 4)    # A right 4
    ]
    
    current_board = initial_board
    solution = []
    
    for car, direction in solution_sequence:
        cars = get_car_positions(current_board)
        if car not in cars:
            print(f"Car {car} not found on board")
            return None
            
        positions = cars[car]
        if can_move(current_board, car, positions, direction):
            new_board = make_move(current_board, car, positions, direction)
            if new_board is None:
                print(f"Invalid move for car {car} direction {direction}")
                return None
            current_board = new_board
            move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
            solution.append(move_str)
        else:
            print(f"Cannot move car {car} in direction {direction}")
            print("Current board state:")
            print_board(current_board)
            return None
    
    # Verify final position of red car (A)
    final_cars = get_car_positions(current_board)
    if 'A' in final_cars:
        red_car_positions = final_cars['A']
        rightmost_position = max(pos[1] for pos in red_car_positions)
        if rightmost_position == 4:  # Check if red car reached the exit
            return solution
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")