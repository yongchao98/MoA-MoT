from collections import deque
import copy

def create_initial_board():
    return [
        ['B', 'B', 'B', '.', 'H', 'I'],
        ['.', '.', 'F', '.', 'H', 'I'],
        ['A', 'A', 'F', '.', '.', 'J'],
        ['E', 'C', 'C', '.', '.', 'J'],
        ['E', '.', 'G', '.', '.', 'J'],
        ['x', '.', 'G', 'D', 'D', 'D']
    ]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, car, row, col, is_vertical, length, delta):
    if is_vertical:
        if delta > 0:  # Moving down
            if row + length - 1 + delta >= 6:
                return False
            return all(board[row + length + i][col] == '.' for i in range(delta))
        else:  # Moving up
            if row + delta < 0:
                return False
            return all(board[row + delta + i][col] == '.' for i in range(-delta))
    else:
        if delta > 0:  # Moving right
            if col + length - 1 + delta >= 6:
                return False
            return all(board[row][col + length + i] == '.' for i in range(delta))
        else:  # Moving left
            if col + delta < 0:
                return False
            return all(board[row][col + delta + i] == '.' for i in range(-delta))

def move_car(board, car, row, col, is_vertical, length, delta):
    new_board = copy.deepcopy(board)
    # Clear old position
    for i in range(length):
        if is_vertical:
            new_board[row + i][col] = '.'
        else:
            new_board[row][col + i] = '.'
    # Set new position
    for i in range(length):
        if is_vertical:
            new_board[row + delta + i][col] = car
        else:
            new_board[row][col + delta + i] = car
    return new_board

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                car = board[i][j]
                if car not in cars:
                    # Find orientation and length
                    is_vertical = (i + 1 < 6 and board[i + 1][j] == car)
                    length = 1
                    if is_vertical:
                        while i + length < 6 and board[i + length][j] == car:
                            length += 1
                    else:
                        while j + length < 6 and board[i][j + length] == car:
                            length += 1
                    cars[car] = (i, j, is_vertical, length)
    return cars

def solve_puzzle():
    initial_board = create_initial_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        if current_board[2][4] == 'A':
            return moves
        
        cars = get_car_info(current_board)
        for car, (row, col, is_vertical, length) in cars.items():
            # Try all possible moves (-3 to +3)
            for delta in range(-3, 4):
                if delta == 0:
                    continue
                if is_valid_move(current_board, car, row, col, is_vertical, length, delta):
                    new_board = move_car(current_board, car, row, col, is_vertical, length, delta)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if delta > 0 else ''}{delta}"]
                        queue.append((new_board, new_moves))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')