from collections import deque
import copy

def create_board():
    return [
        ['B', 'B', 'C', 'C', '.', 'x'],
        ['.', '.', '.', 'K', 'D', 'D'],
        ['x', 'A', 'A', 'K', '.', 'L'],
        ['E', 'E', 'J', 'F', 'F', 'L'],
        ['.', '.', 'J', '.', 'G', 'G'],
        ['.', 'H', 'H', 'H', 'I', 'I']
    ]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
                car = board[i][j]
                if car not in cars:
                    cars[car] = []
                cars[car].append((i, j))
    return cars

def is_valid_move(board, car_pos, direction, distance):
    if len(car_pos) < 2:
        return False
    
    is_horizontal = car_pos[0][0] == car_pos[1][0]
    row, col = car_pos[0]
    
    if is_horizontal:
        if direction not in ['+', '-']:
            return False
        delta = 1 if direction == '+' else -1
        for pos in car_pos:
            for d in range(1, abs(distance) + 1):
                new_col = pos[1] + d * delta
                if new_col < 0 or new_col >= 6 or board[pos[0]][new_col] not in ['.']:
                    return False
    else:
        if direction not in ['+', '-']:
            return False
        delta = 1 if direction == '+' else -1
        for pos in car_pos:
            for d in range(1, abs(distance) + 1):
                new_row = pos[0] + d * delta
                if new_row < 0 or new_row >= 6 or board[new_row][pos[1]] not in ['.']:
                    return False
    return True

def make_move(board, car, car_pos, direction, distance):
    new_board = [row[:] for row in board]
    is_horizontal = car_pos[0][0] == car_pos[1][0]
    delta = 1 if direction == '+' else -1
    
    # Clear current positions
    for row, col in car_pos:
        new_board[row][col] = '.'
    
    # Set new positions
    for row, col in car_pos:
        if is_horizontal:
            new_board[row][col + delta * distance] = car
        else:
            new_board[row + delta * distance][col] = car
            
    return new_board

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        cars = get_car_positions(current_board)
        if 'A' in cars and any(col == 4 for row, col in cars['A']):
            return moves
        
        # Try all possible moves for each car
        for car, positions in cars.items():
            for direction in ['+', '-']:
                for distance in range(1, 6):
                    if is_valid_move(current_board, positions, direction, distance):
                        new_board = make_move(current_board, car, positions, direction, distance)
                        board_str = board_to_string(new_board)
                        
                        if board_str not in visited:
                            visited.add(board_str)
                            new_moves = moves + [f"{car}{direction}{distance}"]
                            queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")