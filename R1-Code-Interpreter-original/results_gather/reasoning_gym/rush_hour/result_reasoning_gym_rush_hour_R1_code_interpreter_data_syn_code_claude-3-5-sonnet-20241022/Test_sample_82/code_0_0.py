from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def is_valid_move(board, car, direction, steps):
    rows, cols = len(board), len(board[0])
    car_positions = [(i, j) for i in range(rows) for j in range(cols) if board[i][j] == car]
    
    # Determine if car is horizontal or vertical
    is_horizontal = car_positions[0][0] == car_positions[1][0]
    
    # Check if move is valid
    for pos in car_positions:
        row, col = pos
        if is_horizontal:
            new_col = col + steps
            if new_col < 0 or new_col >= cols:
                return False
            if steps > 0:
                for i in range(1, steps + 1):
                    if col + i >= cols or board[row][col + i] not in ['.', car]:
                        return False
            else:
                for i in range(1, -steps + 1):
                    if col - i < 0 or board[row][col - i] not in ['.', car]:
                        return False
        else:
            new_row = row + steps
            if new_row < 0 or new_row >= rows:
                return False
            if steps > 0:
                for i in range(1, steps + 1):
                    if row + i >= rows or board[row + i][col] not in ['.', car]:
                        return False
            else:
                for i in range(1, -steps + 1):
                    if row - i < 0 or board[row - i][col] not in ['.', car]:
                        return False
    return True

def move_car(board, car, steps):
    new_board = [list(row) for row in board]
    rows, cols = len(board), len(board[0])
    car_positions = [(i, j) for i in range(rows) for j in range(cols) if board[i][j] == car]
    
    # Clear car from current position
    for row, col in car_positions:
        new_board[row][col] = '.'
    
    # Place car in new position
    is_horizontal = car_positions[0][0] == car_positions[1][0]
    if is_horizontal:
        row = car_positions[0][0]
        start_col = min(pos[1] for pos in car_positions)
        for i in range(len(car_positions)):
            new_board[row][start_col + steps + i] = car
    else:
        col = car_positions[0][1]
        start_row = min(pos[0] for pos in car_positions)
        for i in range(len(car_positions)):
            new_board[start_row + steps + i][col] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        'F..HBB',
        'FCCH.I',
        'AAG..I',
        '..GDDI',
        '..xx..',
        '...EEE'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) has reached the exit
        for i in range(len(current_board)):
            if 'AA' in current_board[i] and current_board[i].rindex('A') == len(current_board[i])-1:
                return moves
        
        # Try all possible moves for all cars
        cars = set(''.join(current_board)) - {'.', 'x'}
        for car in cars:
            for steps in range(-5, 6):
                if steps == 0:
                    continue
                if is_valid_move(current_board, car, 'horizontal', steps):
                    new_board = move_car(current_board, car, steps)
                    board_tuple = tuple(new_board)
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        queue.append((new_board, moves + [f"{car}{'+' if steps > 0 else ''}{steps}"]))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")