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
    
    for car in cars:
        pos = cars[car]['positions']
        pos.sort()
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[-1][0] else 'V'
    return cars

def get_blocking_cars(board):
    # Find red car position
    red_row = None
    red_col = None
    for i in range(6):
        for j in range(6):
            if board[i][j] == 'A':
                red_row = i
                red_col = j + 1
                break
        if red_row is not None:
            break
    
    # Get cars blocking the path
    blocking = set()
    for j in range(red_col + 1, 6):
        if board[red_row][j] != '.' and board[red_row][j] != 'x':
            blocking.add(board[red_row][j])
    return blocking

def is_valid_move(board, positions, orientation, distance):
    if orientation == 'H':
        row = positions[0][0]
        if distance > 0:  # right
            return all(j + distance < 6 and board[row][j + distance] == '.' 
                      for j in range(positions[-1][1] + 1, positions[-1][1] + distance + 1))
        else:  # left
            return all(j + distance >= 0 and board[row][j + distance] == '.' 
                      for j in range(positions[0][1] + distance, positions[0][1]))
    else:  # vertical
        col = positions[0][1]
        if distance > 0:  # down
            return all(i + distance < 6 and board[i + distance][col] == '.' 
                      for i in range(positions[-1][0] + 1, positions[-1][0] + distance + 1))
        else:  # up
            return all(i + distance >= 0 and board[i + distance][col] == '.' 
                      for i in range(positions[0][0] + distance, positions[0][0]))

def make_move(board, car, positions, orientation, distance):
    new_board = [list(row) for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    if orientation == 'H':
        row = positions[0][0]
        for j in range(len(positions)):
            new_board[row][positions[j][1] + distance] = car
    else:
        col = positions[0][1]
        for i in range(len(positions)):
            new_board[positions[i][0] + distance][col] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(board):
    initial_state = tuple(board)
    visited = {initial_state}
    queue = deque([(initial_state, [])])
    max_moves = 10  # Limit search depth
    
    while queue:
        current_state, moves = queue.popleft()
        current_board = list(current_state)
        
        if len(moves) > max_moves:
            continue
            
        if 'A' in current_board[2][-2:]:
            return moves
            
        blocking_cars = get_blocking_cars(current_board)
        cars = get_car_info(current_board)
        
        # Prioritize moves of blocking cars and the red car
        priority_cars = ['A'] + list(blocking_cars)
        for car in priority_cars:
            if car not in cars:
                continue
                
            info = cars[car]
            positions = info['positions']
            orientation = info['orientation']
            
            # Try minimal moves first
            for distance in [-2, -1, 1, 2]:
                if is_valid_move(current_board, positions, orientation, distance):
                    new_board = make_move(current_board, car, positions, orientation, distance)
                    new_state = tuple(new_board)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{car}{'+' if distance > 0 else '-'}{abs(distance)}"]
                        queue.append((new_state, new_moves))
        
        # If no solution found with priority cars, try others
        if len(queue) == 0:
            max_moves += 5  # Increase search depth
            queue.append((initial_state, []))
            visited = {initial_state}
    
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