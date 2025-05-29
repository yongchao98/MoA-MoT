from collections import deque

def get_board_string(board):
    return ''.join(''.join(row) for row in board)

def get_blocking_cars(board, red_car_row):
    # Get cars blocking the red car's path to exit
    blocking = []
    right_edge = -1
    # Find the rightmost position of the red car
    for j in range(len(board[0])):
        if board[red_car_row][j] == 'A':
            right_edge = j
    
    # Check all cars to the right of the red car
    for j in range(right_edge + 1, len(board[0])):
        if board[red_car_row][j] not in ['.', 'x'] and board[red_car_row][j] not in blocking:
            blocking.append(board[red_car_row][j])
    
    return blocking

def can_move(board, car, positions, direction):
    is_horizontal = positions[0][0] == positions[1][0]
    if is_horizontal:
        row = positions[0][0]
        if direction < 0:  # Try left
            col = min(p[1] for p in positions) + direction
            return col >= 0 and board[row][col] == '.'
        else:  # Try right
            col = max(p[1] for p in positions) + direction
            return col < len(board[0]) and board[row][col] == '.'
    else:
        col = positions[0][1]
        if direction < 0:  # Try up
            row = min(p[0] for p in positions) + direction
            return row >= 0 and board[row][col] == '.'
        else:  # Try down
            row = max(p[0] for p in positions) + direction
            return row < len(board) and board[row][col] == '.'

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    for pos in positions:
        if is_horizontal:
            new_board[pos[0]][pos[1] + direction] = car
        else:
            new_board[pos[0] + direction][pos[1]] = car
            
    return [''.join(row) for row in new_board]

def solve_rush_hour():
    initial_board = [
        "BBBKCC",
        "DDJK.L",
        "I.JAAL",
        "IEE.xM",
        "FF...M",
        "GG.xHH"
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(get_board_string(initial_board))
    
    # Find red car row
    red_car_row = next(i for i, row in enumerate(initial_board) if 'A' in row)
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved
        if 'A' in cars and max(p[1] for p in cars['A']) == len(current_board[0]) - 1:
            return moves
        
        # Get blocking cars
        blocking = get_blocking_cars(current_board, red_car_row)
        
        # Prioritize moves
        for car in ['A'] + blocking + [c for c in cars if c not in blocking and c != 'A']:
            positions = cars[car]
            is_horizontal = positions[0][0] == positions[1][0]
            
            for direction in [-1, 1]:
                if can_move(current_board, car, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_string = get_board_string(new_board)
                    
                    if board_string not in visited:
                        visited.add(board_string)
                        new_moves = moves + [(car, direction)]
                        
                        # Prioritize promising moves
                        if car == 'A' and direction > 0:
                            queue.appendleft((new_board, new_moves))
                        elif car in blocking:
                            queue.append((new_board, new_moves))
                        else:
                            queue.append((new_board, new_moves))
    
    return None

# Solve and format the output
solution = solve_rush_hour()
if solution:
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}1")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("No solution found")