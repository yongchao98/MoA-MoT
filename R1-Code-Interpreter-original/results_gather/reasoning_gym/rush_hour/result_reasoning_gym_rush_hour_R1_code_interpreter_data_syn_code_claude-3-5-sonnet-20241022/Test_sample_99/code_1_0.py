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
    return {k: sorted(v) for k, v in cars.items()}

def is_valid_move(board, car, car_pos, direction, distance):
    is_horizontal = car_pos[0][0] == car_pos[1][0]
    row, col = car_pos[0]
    
    if is_horizontal:
        if direction not in ['+', '-']:
            return False
        delta = 1 if direction == '+' else -1
        for pos in car_pos:
            for d in range(1, distance + 1):
                new_col = pos[1] + (delta * d)
                if new_col < 0 or new_col >= 6:
                    return False
                if board[pos[0]][new_col] != '.' and board[pos[0]][new_col] != 'x':
                    return False
                if board[pos[0]][new_col] == 'x' and car != 'A':
                    return False
    else:
        if direction not in ['+', '-']:
            return False
        delta = 1 if direction == '+' else -1
        for pos in car_pos:
            for d in range(1, distance + 1):
                new_row = pos[0] + (delta * d)
                if new_row < 0 or new_row >= 6:
                    return False
                if board[new_row][pos[1]] != '.':
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
            new_col = col + (delta * distance)
            new_board[row][new_col] = car
        else:
            new_row = row + (delta * distance)
            new_board[new_row][col] = car
            
    return new_board

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (AA is at the exit)
        if 'A' in cars:
            a_positions = cars['A']
            if any(col == 4 for row, col in a_positions):
                return moves
        
        # Try all possible moves for each car
        for car, positions in cars.items():
            for direction in ['+', '-']:
                for distance in range(1, 4):  # Limit distance to 3 spaces
                    if is_valid_move(current_board, car, positions, direction, distance):
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