from collections import deque

def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            car = board[i][j]
            if car != '.' and car != 'x':
                if car not in cars:
                    cars[car] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[car]['positions'].append((i, j))
                    # Determine orientation
                    if len(cars[car]['positions']) == 2:
                        cars[car]['orientation'] = 'H' if cars[car]['positions'][0][0] == cars[car]['positions'][1][0] else 'V'
    return cars

def get_blocking_cars(cars):
    # Find cars blocking the red car's path
    red_row = cars['A']['positions'][0][0]
    red_right = max(pos[1] for pos in cars['A']['positions'])
    blocking = set()
    
    for car, info in cars.items():
        if car != 'A':
            for pos in info['positions']:
                if pos[0] == red_row and pos[1] > red_right:
                    blocking.add(car)
    return blocking

def can_move(board, car_info, car, direction):
    positions = car_info[car]['positions']
    orientation = car_info[car]['orientation']
    
    if orientation == 'H':
        row = positions[0][0]
        if direction > 0:  # Try right
            rightmost = max(pos[1] for pos in positions)
            return rightmost + 1 < len(board[0]) and board[row][rightmost + 1] == '.'
        else:  # Try left
            leftmost = min(pos[1] for pos in positions)
            return leftmost > 0 and board[row][leftmost - 1] == '.'
    else:  # Vertical
        col = positions[0][1]
        if direction > 0:  # Try down
            bottommost = max(pos[0] for pos in positions)
            return bottommost + 1 < len(board) and board[bottommost + 1][col] == '.'
        else:  # Try up
            topmost = min(pos[0] for pos in positions)
            return topmost > 0 and board[topmost - 1][col] == '.'

def move_car(board, car, direction):
    new_board = [list(row) for row in board]
    car_positions = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == car]
    is_horizontal = car_positions[0][0] == car_positions[1][0]
    
    # Clear current positions
    for i, j in car_positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in car_positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "BBBCC.",
        "..JKDD",
        "AAJKLM",
        "I.EELM",
        "IFF..N",
        "GGHHHN"
    ]
    
    visited = set(tuple(initial_board))
    queue = deque([(initial_board, [])])
    
    while queue:
        current_board, moves = queue.popleft()
        car_info = get_car_info(current_board)
        blocking_cars = get_blocking_cars(car_info)
        
        # Check if solved
        red_car = car_info['A']
        if max(pos[1] for pos in red_car['positions']) == len(current_board[0]) - 1:
            return moves
        
        # Prioritize moves that help clear the path
        for car in blocking_cars:
            for direction in [-1, 1]:
                if can_move(current_board, car_info, car, direction):
                    new_board = move_car(current_board, car, direction)
                    board_tuple = tuple(new_board)
                    if board_tuple not in visited:
                        visited.add(board_tuple)
                        move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                        queue.append((new_board, moves + [move_str]))
        
        # Try other moves
        for car in car_info:
            if car not in blocking_cars:
                for direction in [-1, 1]:
                    if can_move(current_board, car_info, car, direction):
                        new_board = move_car(current_board, car, direction)
                        board_tuple = tuple(new_board)
                        if board_tuple not in visited:
                            visited.add(board_tuple)
                            move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                            queue.append((new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")