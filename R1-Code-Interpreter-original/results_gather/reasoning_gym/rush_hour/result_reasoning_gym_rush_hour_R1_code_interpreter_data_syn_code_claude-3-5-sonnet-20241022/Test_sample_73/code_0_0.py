from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def is_valid_move(board, car, direction, steps):
    rows, cols = len(board), len(board[0])
    # Find car positions
    car_positions = []
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == car:
                car_positions.append((i, j))
    
    # Determine orientation (horizontal or vertical)
    is_horizontal = car_positions[0][0] == car_positions[-1][0]
    
    # Check if move is valid
    if is_horizontal:
        row = car_positions[0][0]
        min_col = min(p[1] for p in car_positions)
        max_col = max(p[1] for p in car_positions)
        if direction > 0:  # Moving right
            new_max = max_col + direction
            if new_max >= cols:
                return False
            return all(board[row][max_col+1:new_max+1] == '.' * direction)
        else:  # Moving left
            new_min = min_col + direction
            if new_min < 0:
                return False
            return all(board[row][new_min:min_col] == '.' * (-direction))
    else:
        col = car_positions[0][1]
        min_row = min(p[0] for p in car_positions)
        max_row = max(p[0] for p in car_positions)
        if direction > 0:  # Moving down
            new_max = max_row + direction
            if new_max >= rows:
                return False
            return all(board[max_row+1:new_max+1][col] == '.' * direction)
        else:  # Moving up
            new_min = min_row + direction
            if new_min < 0:
                return False
            return all(board[new_min:min_row][col] == '.' * (-direction))

def move_car(board, car, direction):
    new_board = [list(row) for row in board]
    rows, cols = len(board), len(board[0])
    # Find car positions
    car_positions = []
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == car:
                car_positions.append((i, j))
                new_board[i][j] = '.'
    
    # Move car
    for i, j in car_positions:
        if car_positions[0][0] == car_positions[-1][0]:  # horizontal
            new_board[i][j + direction] = car
        else:  # vertical
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
        
        # Check if solved (AA is at the exit)
        if 'AA' in current_board[2] and current_board[2].index('A') >= 4:
            return moves
        
        # Try all possible moves for each car
        cars = set(''.join(current_board))
        cars.discard('.')
        cars.discard('x')
        
        for car in cars:
            # Try moving left/up (-1) and right/down (+1)
            for direction in [-1, 1]:
                # Try moving 1 to 3 spaces
                for steps in range(1, 4):
                    if is_valid_move(current_board, car, direction * steps):
                        new_board = move_car(current_board, car, direction * steps)
                        board_tuple = tuple(new_board)
                        
                        if board_tuple not in seen:
                            seen.add(board_tuple)
                            new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{steps}"]
                            queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")