from collections import deque

def get_board_state(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, car, direction, distance):
    rows, cols = len(board), len(board[0])
    car_positions = [(r, c) for r in range(rows) for c in range(cols) if board[r][c] == car]
    
    # Determine orientation (horizontal or vertical)
    is_horizontal = car_positions[0][0] == car_positions[-1][0]
    
    for d in range(1, abs(distance) + 1):
        for r, c in car_positions:
            new_r = r + (0 if is_horizontal else (d if distance > 0 else -d))
            new_c = c + (d if is_horizontal and distance > 0 else -d if is_horizontal else 0)
            
            if not (0 <= new_r < rows and 0 <= new_c < cols):
                return False
            if board[new_r][new_c] not in ['.', car]:
                return False
    return True

def make_move(board, car, direction, distance):
    new_board = [list(row) for row in board]
    rows, cols = len(board), len(board[0])
    car_positions = [(r, c) for r in range(rows) for c in range(cols) if board[r][c] == car]
    
    # Clear current positions
    for r, c in car_positions:
        new_board[r][c] = '.'
    
    # Determine orientation
    is_horizontal = car_positions[0][0] == car_positions[-1][0]
    
    # Place car in new positions
    for r, c in car_positions:
        new_r = r + (0 if is_horizontal else distance)
        new_c = c + (distance if is_horizontal else 0)
        new_board[new_r][new_c] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    visited = {get_board_state(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) is at exit
        for row in current_board:
            if row.endswith('AA'):
                return moves
        
        # Try all possible moves for each car
        cars = set(c for row in current_board for c in row if c.isalpha())
        for car in cars:
            for distance in [-3, -2, -1, 1, 2, 3]:
                if is_valid_move(current_board, car, 'horizontal' if car == car else 'vertical', distance):
                    new_board = make_move(current_board, car, 'horizontal' if car == car else 'vertical', distance)
                    new_state = get_board_state(new_board)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{car}{'+' if distance > 0 else ''}{distance}"]
                        queue.append((new_board, new_moves))
    
    return None

# Initial board
initial_board = [
    '.EBBB.',
    '.E..F.',
    'AA..F.',
    '..CCC.',
    '..DDD.',
    '......'
]

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print('No solution found')