from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n')]

def get_blocking_cars(board, car_pos):
    """Identify cars directly blocking the path to exit for car A"""
    blocking = set()
    row = car_pos[0][0]  # Row of car A
    rightmost = max(p[1] for p in car_pos)
    
    # Check all positions to the right of car A
    for col in range(rightmost + 1, 6):
        if board[row][col] not in '.x':
            blocking.add(board[row][col])
    return blocking

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'positions': [(i, j)]}
                else:
                    cars[board[i][j]]['positions'].append((i, j))
    
    # Determine orientation
    for car in cars:
        cars[car]['horizontal'] = cars[car]['positions'][0][0] == cars[car]['positions'][1][0]
    return cars

def can_move(board, car_info, direction):
    positions = car_info['positions']
    if car_info['horizontal']:
        row = positions[0][0]
        if direction < 0:  # Left
            col = min(p[1] for p in positions)
            return col > 0 and board[row][col-1] == '.'
        else:  # Right
            col = max(p[1] for p in positions)
            return col < 5 and board[row][col+1] == '.'
    else:
        col = positions[0][1]
        if direction < 0:  # Up
            row = min(p[0] for p in positions)
            return row > 0 and board[row-1][col] == '.'
        else:  # Down
            row = max(p[0] for p in positions)
            return row < 5 and board[row+1][col] == '.'

def make_move(board, car, car_info, direction):
    new_board = [row[:] for row in board]
    positions = car_info['positions']
    horizontal = car_info['horizontal']
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    new_positions = []
    for i, j in positions:
        new_i = i + (0 if horizontal else direction)
        new_j = j + (direction if horizontal else 0)
        new_board[new_i][new_j] = car
        new_positions.append((new_i, new_j))
    
    return new_board, new_positions

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    cars = get_car_positions(board)
    
    queue = deque([(board, cars, [])])
    visited = {tuple(map(tuple, board))}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if solved
        red_car = current_cars['A']
        if max(pos[1] for pos in red_car['positions']) == 5:
            return moves
        
        # Get blocking cars
        blocking = get_blocking_cars(current_board, red_car['positions'])
        
        # Prioritize moving blocking cars and car A
        priority_cars = list(blocking) + ['A']
        all_cars = priority_cars + [c for c in current_cars if c not in priority_cars]
        
        for car in all_cars:
            car_info = current_cars[car]
            for direction in [-1, 1]:
                if can_move(current_board, car_info, direction):
                    new_board, new_positions = make_move(current_board, car, car_info, direction)
                    board_tuple = tuple(map(tuple, new_board))
                    
                    if board_tuple not in visited:
                        visited.add(board_tuple)
                        new_cars = dict(current_cars)
                        new_cars[car] = {
                            'positions': new_positions,
                            'horizontal': car_info['horizontal']
                        }
                        new_moves = moves + [(car, direction)]
                        queue.append((new_board, new_cars, new_moves))
    
    return None

# Initial board
initial_board = """..xBBK
GCCJ.K
GAAJ..
DDIJ..
.HIEE.
.HFFx."""

# Find and format solution
solution = solve_puzzle(initial_board)
if solution:
    move_str = ''
    for car, direction in solution:
        move_str += f'{car}{"+"+str(direction) if direction > 0 else str(direction)} '
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")