from collections import deque

def get_car_info(board):
    # Returns dict of car positions and their orientations (h: horizontal, v: vertical)
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                car = board[i][j]
                if car not in cars:
                    cars[car] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[car]['positions'].append((i, j))
                    # Determine orientation
                    if len(cars[car]['positions']) == 2:
                        if cars[car]['positions'][0][0] == cars[car]['positions'][1][0]:
                            cars[car]['orientation'] = 'h'  # horizontal
                        else:
                            cars[car]['orientation'] = 'v'  # vertical
    return cars

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_blocking_cars(board, cars):
    # Get cars blocking the red car's path
    aa_info = cars['A']
    aa_row = aa_info['positions'][0][0]
    aa_rightmost = max(pos[1] for pos in aa_info['positions'])
    
    blocking = []
    for col in range(aa_rightmost + 1, len(board[0])):
        if board[aa_row][col] not in ['.', 'x']:
            blocking.append(board[aa_row][col])
    return blocking

def can_move(board, car_info, direction):
    positions = car_info['positions']
    orientation = car_info['orientation']
    
    if orientation == 'h':
        row = positions[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_info, direction):
    new_board = [row[:] for row in board]
    positions = car_info['positions']
    orientation = car_info['orientation']
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if orientation == 'h':
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = car
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = car
            
    return new_board

def solve_puzzle():
    initial_board = [
        ['x','.','G','B','B','J'],
        ['C','C','G','H','I','J'],
        ['F','A','A','H','I','K'],
        ['F','D','D','.','I','K'],
        ['E','E','.','x','.','.'],
        ['.','.','.','.','.','.']
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_info(current_board)
        
        # Check if solved
        if 'A' in cars:
            aa_pos = cars['A']['positions']
            if max(pos[1] for pos in aa_pos) == len(current_board[0])-1:
                return moves
        
        # Get blocking cars and prioritize moving them
        blocking_cars = get_blocking_cars(current_board, cars)
        priority_cars = blocking_cars + [car for car in cars if car not in blocking_cars]
        
        for car in priority_cars:
            car_info = cars[car]
            for direction in [-1, 1]:
                if can_move(current_board, car_info, direction):
                    new_board = move_car(current_board, car, car_info, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")