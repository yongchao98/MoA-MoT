from collections import deque

def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_valid_move(board, car_positions, direction, is_horizontal):
    rows, cols = len(board), len(board[0])
    if is_horizontal:
        row = car_positions[0][0]
        min_col = min(j for _, j in car_positions)
        max_col = max(j for _, j in car_positions)
        if direction > 0:  # Moving right
            new_max = max_col + direction
            if new_max >= cols: return False
            return all(board[row][j] == '.' for j in range(max_col + 1, new_max + 1))
        else:  # Moving left
            new_min = min_col + direction
            if new_min < 0: return False
            return all(board[row][j] == '.' for j in range(new_min, min_col))
    else:
        col = car_positions[0][1]
        min_row = min(i for i, _ in car_positions)
        max_row = max(i for i, _ in car_positions)
        if direction > 0:  # Moving down
            new_max = max_row + direction
            if new_max >= rows: return False
            return all(board[i][col] == '.' for i in range(max_row + 1, new_max + 1))
        else:  # Moving up
            new_min = min_row + direction
            if new_min < 0: return False
            return all(board[i][col] == '.' for i in range(new_min, min_row))

def move_car(board, car, car_positions, direction, is_horizontal):
    new_board = [list(row) for row in board]
    # Clear old positions
    for i, j in car_positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in car_positions:
        if is_horizontal:
            new_board[i][j + direction] = car
        else:
            new_board[i + direction][j] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        'HIBBCC',
        'HIJ...',
        'AAJ..L',
        'DDJK.L',
        'xEEK.L',
        '.FFGG.'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Get current state of all cars
        cars = get_car_info(current_board)
        
        # Check if solved (AA is at the exit)
        if 'A' in cars:
            a_positions = cars['A']
            if max(j for _, j in a_positions) == 5:  # Red car reached the exit
                return moves
        
        # Process each car
        for car, positions in cars.items():
            is_horizontal = positions[0][0] == positions[-1][0]
            
            # Determine possible move distances
            if is_horizontal:
                max_distance = 4 if car == 'A' else 3
                distances = list(range(-max_distance, max_distance + 1))
            else:
                distances = [-2, -1, 1, 2]
            
            # Try each possible move
            for d in distances:
                if d != 0 and is_valid_move(current_board, positions, d, is_horizontal):
                    new_board = move_car(current_board, car, positions, d, is_horizontal)
                    board_tuple = tuple(new_board)
                    
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        new_moves = moves + [f"{car}{'+' if d > 0 else '-'}{abs(d)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")