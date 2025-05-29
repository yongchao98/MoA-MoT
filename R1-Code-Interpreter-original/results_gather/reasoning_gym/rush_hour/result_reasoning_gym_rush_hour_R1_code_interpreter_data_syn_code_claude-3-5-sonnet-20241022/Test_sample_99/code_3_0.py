from collections import deque

def create_board():
    return [
        list('BBCC.x'),
        list('...KDD'),
        list('xAAK.L'),
        list('EEJFFL'),
        list('..J.GG'),
        list('.HHHII')
    ]

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

def get_car_orientation(positions):
    if len(positions) < 2:
        return None
    # If all positions have same row, it's horizontal
    if all(pos[0] == positions[0][0] for pos in positions):
        return 'H'
    # If all positions have same column, it's vertical
    if all(pos[1] == positions[0][1] for pos in positions):
        return 'V'
    return None

def can_move(board, car, car_positions, direction):
    orientation = get_car_orientation(car_positions)
    if not orientation:
        return False
    
    if orientation == 'H':
        row = car_positions[0][0]
        if direction == '+':
            max_col = max(col for _, col in car_positions)
            if max_col + 1 >= 6:
                return False
            return (board[row][max_col + 1] == '.' or 
                   (car == 'A' and board[row][max_col + 1] == 'x'))
        else:
            min_col = min(col for _, col in car_positions)
            if min_col - 1 < 0:
                return False
            return board[row][min_col - 1] == '.'
    else:  # Vertical
        col = car_positions[0][1]
        if direction == '+':
            max_row = max(row for row, _ in car_positions)
            if max_row + 1 >= 6:
                return False
            return board[max_row + 1][col] == '.'
        else:
            min_row = min(row for row, _ in car_positions)
            if min_row - 1 < 0:
                return False
            return board[min_row - 1][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [row[:] for row in board]
    orientation = get_car_orientation(positions)
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if orientation == 'H':
        delta = 1 if direction == '+' else -1
        row = positions[0][0]
        for _, col in positions:
            new_board[row][col + delta] = car
    else:
        delta = 1 if direction == '+' else -1
        col = positions[0][1]
        for row, _ in positions:
            new_board[row + delta][col] = car
            
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

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
            if max(col for _, col in a_positions) == 4:
                return moves
        
        # Try moves for each car
        for car, positions in cars.items():
            for direction in ['+', '-']:
                if can_move(current_board, car, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{direction}1"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")