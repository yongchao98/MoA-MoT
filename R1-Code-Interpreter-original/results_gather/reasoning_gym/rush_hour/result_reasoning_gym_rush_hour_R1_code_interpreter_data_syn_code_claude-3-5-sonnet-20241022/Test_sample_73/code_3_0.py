from collections import deque

def is_valid_move(board, car, direction):
    rows, cols = len(board), len(board[0])
    car_positions = []
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == car:
                car_positions.append((i, j))
    
    is_horizontal = car_positions[0][0] == car_positions[-1][0]
    
    if is_horizontal:
        row = car_positions[0][0]
        min_col = min(p[1] for p in car_positions)
        max_col = max(p[1] for p in car_positions)
        if direction > 0:
            new_max = max_col + direction
            if new_max >= cols: return False
            return all(board[row][j] == '.' for j in range(max_col + 1, new_max + 1))
        else:
            new_min = min_col + direction
            if new_min < 0: return False
            return all(board[row][j] == '.' for j in range(new_min, min_col))
    else:
        col = car_positions[0][1]
        min_row = min(p[0] for p in car_positions)
        max_row = max(p[0] for p in car_positions)
        if direction > 0:
            new_max = max_row + direction
            if new_max >= rows: return False
            return all(board[i][col] == '.' for i in range(max_row + 1, new_max + 1))
        else:
            new_min = min_row + direction
            if new_min < 0: return False
            return all(board[i][col] == '.' for i in range(new_min, min_row))

def move_car(board, car, direction):
    new_board = [list(row) for row in board]
    car_positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == car:
                car_positions.append((i, j))
                new_board[i][j] = '.'
    
    is_horizontal = car_positions[0][0] == car_positions[-1][0]
    for i, j in car_positions:
        if is_horizontal:
            new_board[i][j + direction] = car
        else:
            new_board[i + direction][j] = car
    
    return [''.join(row) for row in new_board]

def get_car_orientation(board, car):
    positions = [(i, j) for i in range(len(board)) for j in range(len(board[0])) if board[i][j] == car]
    return 'horizontal' if positions[0][0] == positions[-1][0] else 'vertical'

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
    
    # Priority order for cars to move
    blocking_cars = ['J', 'K', 'L']  # Cars directly blocking AA
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        if 'AA' in current_board[2] and current_board[2].index('A') >= 4:
            return moves
        
        # Try moving blocking cars
        for car in blocking_cars:
            orientation = get_car_orientation(current_board, car)
            directions = [-1, 1] if orientation == 'vertical' else [-3, -2, -1, 1, 2, 3]
            
            for direction in directions:
                if is_valid_move(current_board, car, direction):
                    new_board = move_car(current_board, car, direction)
                    board_tuple = tuple(new_board)
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
        
        # Try moving the red car
        if is_valid_move(current_board, 'A', 4):  # Try moving all the way to exit
            new_board = move_car(current_board, 'A', 4)
            board_tuple = tuple(new_board)
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + ['A+4']
                queue.append((new_board, new_moves))
        elif is_valid_move(current_board, 'A', 1):  # Try moving one step
            new_board = move_car(current_board, 'A', 1)
            board_tuple = tuple(new_board)
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + ['A+1']
                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")