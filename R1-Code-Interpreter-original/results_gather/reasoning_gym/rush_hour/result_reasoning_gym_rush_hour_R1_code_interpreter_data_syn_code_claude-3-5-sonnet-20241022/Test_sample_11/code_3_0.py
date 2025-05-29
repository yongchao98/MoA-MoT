from collections import deque

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i,j)]
                else:
                    cars[board[i][j]].append((i,j))
    return cars

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        # Check if horizontal (same row)
        if positions[0][0] == positions[-1][0]:
            row = positions[0][0]
            min_col = min(p[1] for p in positions)
            max_col = max(p[1] for p in positions)
            # Try move left
            if min_col > 0 and board[row][min_col-1] == '.':
                moves.append((car, -1))
            # Try move right
            if max_col < 5 and board[row][max_col+1] == '.':
                moves.append((car, 1))
        else:  # Vertical car
            col = positions[0][1]
            min_row = min(p[0] for p in positions)
            max_row = max(p[0] for p in positions)
            # Try move up
            if min_row > 0 and board[min_row-1][col] == '.':
                moves.append((car, -1))
            # Try move down
            if max_row < 5 and board[max_row+1][col] == '.':
                moves.append((car, 1))
    return moves

def make_move(board, car, direction):
    board = [list(row) for row in board]
    cars = get_car_positions(board)
    positions = sorted(cars[car])
    
    # Determine orientation
    is_horizontal = positions[0][0] == positions[-1][0]
    
    # Remove car from current position
    for i, j in positions:
        board[i][j] = '.'
    
    # Add car in new position
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        board[new_i][new_j] = car
    
    return [''.join(row) for row in board]

def can_reach_exit(board):
    cars = get_car_positions(board)
    red_car = cars['A']
    # Check if path to exit is clear
    row = red_car[0][0]
    max_col = max(p[1] for p in red_car)
    for j in range(max_col + 1, 6):
        if board[row][j] != '.':
            return False
    return True

def solve_puzzle():
    initial_board = [
        'EBBCCC',
        'E....H',
        'F.xAAH',
        'F.G...',
        '..GDDD',
        '......'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if we can reach the exit
        if can_reach_exit(current_board):
            cars = get_car_positions(current_board)
            red_car = cars['A']
            right_pos = max(p[1] for p in red_car)
            # Add final move to exit
            return moves + [('A', 6 - right_pos)]
        
        cars = get_car_positions(current_board)
        for car, direction in get_valid_moves(current_board, cars):
            new_board = make_move(current_board, car, direction)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}{abs(direction)}")
    print('<<<' + ' '.join(formatted_moves) + '>>>')
else:
    print("No solution found")