from collections import deque

def create_board():
    return [
        ['B', 'B', 'B', 'K', 'C', 'C'],
        ['D', 'D', 'J', 'K', '.', 'L'],
        ['I', '.', 'J', 'A', 'A', 'L'],
        ['I', 'E', 'E', '.', 'x', 'M'],
        ['F', 'F', '.', '.', '.', 'M'],
        ['G', 'G', 'H', 'H', '.', 'x']
    ]

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            car = board[i][j]
            if car not in ['.', 'x']:
                if car not in cars:
                    cars[car] = [(i, j)]
                else:
                    cars[car].append((i, j))
    
    # Determine orientation and sort positions
    for car in cars:
        cars[car].sort()
        is_horizontal = cars[car][0][0] == cars[car][-1][0]
        cars[car] = {
            'positions': cars[car],
            'horizontal': is_horizontal
        }
    return cars

def is_valid_move(board, positions, is_horizontal, delta):
    if is_horizontal:
        row = positions[0][0]
        if delta > 0:  # right
            new_col = positions[-1][1] + 1
            return new_col < 6 and board[row][new_col] == '.'
        else:  # left
            new_col = positions[0][1] - 1
            return new_col >= 0 and board[row][new_col] == '.'
    else:
        col = positions[0][1]
        if delta > 0:  # down
            new_row = positions[-1][0] + 1
            return new_row < 6 and board[new_row][col] == '.'
        else:  # up
            new_row = positions[0][0] - 1
            return new_row >= 0 and board[new_row][col] == '.'

def make_move(board, car, positions, is_horizontal, delta):
    new_board = [row[:] for row in board]
    new_positions = []
    
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in positions:
        new_i = i + (0 if is_horizontal else delta)
        new_j = j + (delta if is_horizontal else 0)
        new_board[new_i][new_j] = car
        new_positions.append((new_i, new_j))
    
    return new_board, new_positions

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    cars = get_car_positions(initial_board)
    
    queue = deque([(initial_board, cars, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if red car (AA) is at exit
        if current_cars['A']['positions'][-1][1] == 4:
            return moves
        
        # Try all possible moves for each car
        for car, info in current_cars.items():
            for delta in [-1, 1]:
                if is_valid_move(current_board, info['positions'], info['horizontal'], delta):
                    new_board, new_positions = make_move(
                        current_board, car, info['positions'], 
                        info['horizontal'], delta
                    )
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        
                        # Update car positions
                        new_cars = dict(current_cars)
                        new_cars[car] = {
                            'positions': new_positions,
                            'horizontal': info['horizontal']
                        }
                        
                        move = f"{car}{'+' if delta > 0 else '-'}1"
                        queue.append((new_board, new_cars, moves + [move]))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print('<<<No solution found>>>')