from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, car_pos, direction):
    if is_horizontal(car_pos):
        if direction > 0:  # Move right
            new_pos = car_pos[-1][1] + 1
            return new_pos < 6 and board[car_pos[0][0]][new_pos] == '.'
        else:  # Move left
            new_pos = car_pos[0][1] - 1
            return new_pos >= 0 and board[car_pos[0][0]][new_pos] == '.'
    else:
        if direction > 0:  # Move down
            new_pos = car_pos[-1][0] + 1
            return new_pos < 6 and board[new_pos][car_pos[0][1]] == '.'
        else:  # Move up
            new_pos = car_pos[0][0] - 1
            return new_pos >= 0 and board[new_pos][car_pos[0][1]] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [list(row) for row in board]
    # Clear old positions
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if is_horizontal(car_pos):
        for i, pos in enumerate(car_pos):
            new_board[pos[0]][pos[1] + direction] = car
    else:
        for i, pos in enumerate(car_pos):
            new_board[pos[0] + direction][pos[1]] = car
            
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        'IBBCC.',
        'IDDEEL',
        'IAA.KL',
        'FFF.K.',
        '.JGGHH',
        '.J....'
    ]
    
    target_x, target_y = 2, 5  # Position where AA needs to reach
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved
        if any(pos == (2, 5) for pos in cars['A']):
            return moves
        
        # Try moving each car
        for car, positions in cars.items():
            # Try both directions
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_tuple = tuple(new_board)
                    
                    if board_tuple not in visited:
                        visited.add(board_tuple)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")