def get_car_positions(board):
    cars = {}
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell not in ['.', 'x']:
                if cell not in cars:
                    cars[cell] = []
                cars[cell].append((i, j))
    return cars

def is_solved(board, red_car_pos):
    return max(x[1] for x in red_car_pos) == len(board[0]) - 1

def get_blocking_path(board, red_car_row, red_car_right):
    blocking = []
    for j in range(red_car_right + 1, len(board[0])):
        if board[red_car_row][j] not in ['.', 'x']:
            blocking.append((board[red_car_row][j], j))
    return blocking

def can_move_vertical(board, car_pos, direction):
    if direction < 0:  # up
        return car_pos[0][0] + direction >= 0 and board[car_pos[0][0] + direction][car_pos[0][1]] == '.'
    else:  # down
        return car_pos[-1][0] + direction < len(board) and board[car_pos[-1][0] + direction][car_pos[-1][1]] == '.'

def can_move_horizontal(board, car_pos, direction):
    if direction < 0:  # left
        return car_pos[0][1] + direction >= 0 and board[car_pos[0][0]][car_pos[0][1] + direction] == '.'
    else:  # right
        return car_pos[-1][1] + direction < len(board[0]) and board[car_pos[0][0]][car_pos[-1][1] + direction] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = car
    
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
    
    # Find red car position
    cars = get_car_positions(initial_board)
    red_car_pos = cars['A']
    red_car_row = red_car_pos[0][0]
    red_car_right = max(x[1] for x in red_car_pos)
    
    visited = set()
    states = [(initial_board, [])]
    visited.add(''.join(initial_board))
    
    while states:
        current_board, moves = states.pop(0)
        cars = get_car_positions(current_board)
        
        # If solved, return moves
        if 'A' in cars and is_solved(current_board, cars['A']):
            return moves
        
        # Get blocking cars in path
        blocking = get_blocking_path(current_board, red_car_row, cars['A'][-1][1])
        
        # Try moving red car first if possible
        if 'A' in cars:
            if can_move_horizontal(current_board, cars['A'], 1):
                new_board = move_car(current_board, 'A', cars['A'], 1)
                board_str = ''.join(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    states.insert(0, (new_board, moves + [('A', 1)]))
                    continue
        
        # Try moving blocking cars
        for car, _ in blocking:
            positions = cars[car]
            is_horizontal = positions[0][0] == positions[1][0]
            
            for direction in [-1, 1]:
                if (is_horizontal and can_move_horizontal(current_board, positions, direction)) or \
                   (not is_horizontal and can_move_vertical(current_board, positions, direction)):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = ''.join(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        states.append((new_board, moves + [(car, direction)]))
        
        # Try moving other cars that might help
        for car, positions in cars.items():
            if car != 'A' and not any(car == b[0] for b in blocking):
                is_horizontal = positions[0][0] == positions[1][0]
                for direction in [-1, 1]:
                    if (is_horizontal and can_move_horizontal(current_board, positions, direction)) or \
                       (not is_horizontal and can_move_vertical(current_board, positions, direction)):
                        new_board = move_car(current_board, car, positions, direction)
                        board_str = ''.join(new_board)
                        if board_str not in visited:
                            visited.add(board_str)
                            states.append((new_board, moves + [(car, direction)]))
    
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