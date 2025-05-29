from collections import deque

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                car = board[i][j]
                if car not in cars:
                    cars[car] = {'positions': [(i, j)], 'length': 1}
                else:
                    cars[car]['positions'].append((i, j))
                    cars[car]['length'] += 1
    
    # Determine orientation and sort positions
    for car in cars:
        pos = cars[car]['positions']
        pos.sort()  # Sort positions for consistent representation
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[-1][0] else 'V'
    
    return cars

def make_move(board, car, car_info, distance):
    new_board = [list(row) for row in board]
    positions = car_info['positions']
    orientation = car_info['orientation']
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    if orientation == 'H':
        row = positions[0][0]
        for col in range(positions[0][1] + distance, positions[0][1] + distance + car_info['length']):
            new_board[row][col] = car
    else:
        col = positions[0][1]
        for row in range(positions[0][0] + distance, positions[0][0] + distance + car_info['length']):
            new_board[row][col] = car
    
    return [''.join(row) for row in new_board]

def is_valid_move(board, car, car_info, distance):
    positions = car_info['positions']
    orientation = car_info['orientation']
    length = car_info['length']
    
    if orientation == 'H':
        row = positions[0][0]
        start_col = positions[0][1]
        
        if distance > 0:  # moving right
            return (start_col + length - 1 + distance < 6 and 
                   all(board[row][start_col + length + d] == '.' for d in range(distance)))
        else:  # moving left
            return (start_col + distance >= 0 and 
                   all(board[row][start_col + d] == '.' for d in range(distance, 0)))
    else:  # vertical
        col = positions[0][1]
        start_row = positions[0][0]
        
        if distance > 0:  # moving down
            return (start_row + length - 1 + distance < 6 and 
                   all(board[start_row + length + d][col] == '.' for d in range(distance)))
        else:  # moving up
            return (start_row + distance >= 0 and 
                   all(board[start_row + d][col] == '.' for d in range(distance, 0)))

def heuristic(board):
    # Simple heuristic: count number of blocking vehicles
    blocking = 0
    red_car_row = None
    red_car_col = None
    
    # Find red car position
    for i in range(6):
        for j in range(6):
            if board[i][j] == 'A':
                red_car_row = i
                red_car_col = j
                break
        if red_car_row is not None:
            break
    
    # Count blocking vehicles
    for j in range(red_car_col + 2, 6):
        if board[red_car_row][j] != '.' and board[red_car_row][j] != 'x':
            blocking += 1
            
    return blocking

def solve_puzzle(board):
    initial_state = tuple(board)
    visited = {initial_state: 0}  # state -> moves count
    queue = [(heuristic(board), 0, initial_state, [])]  # (heuristic, moves_count, state, moves)
    
    while queue:
        queue.sort()  # Sort by heuristic value
        _, _, current_state, moves = queue.pop(0)
        current_board = list(current_state)
        
        if 'A' in current_board[2][-2:]:
            return moves
        
        cars = get_car_info(current_board)
        
        for car in sorted(cars.keys()):  # Sort cars to prioritize certain moves
            info = cars[car]
            
            for distance in range(-5, 6):
                if distance == 0:
                    continue
                    
                if is_valid_move(current_board, car, info, distance):
                    new_board = make_move(current_board, car, info, distance)
                    new_state = tuple(new_board)
                    new_moves = moves + [f"{car}{'+' if distance > 0 else '-'}{abs(distance)}"]
                    
                    if new_state not in visited or len(new_moves) < visited[new_state]:
                        visited[new_state] = len(new_moves)
                        h_score = heuristic(new_board) + len(new_moves)
                        queue.append((h_score, len(new_moves), new_state, new_moves))
    
    return None

# Initial board
board = [
    "GBBIJK",
    "G..IJK",
    "AAHI..",
    "..HCCC",
    "..xDD.",
    "EEEFF."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")