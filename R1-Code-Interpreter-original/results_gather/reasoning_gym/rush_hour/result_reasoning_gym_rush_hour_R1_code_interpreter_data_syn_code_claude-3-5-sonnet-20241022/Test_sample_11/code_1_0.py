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

def get_car_orientation(positions):
    return 'H' if positions[0][0] == positions[-1][0] else 'V'

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        orientation = get_car_orientation(positions)
        if orientation == 'H':
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            # Try left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Try right
            if right < 5 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:  # Vertical
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            # Try up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Try down
            if bottom < 5 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    return moves

def make_move(board, cars, car, direction):
    board = [list(row) for row in board]
    positions = cars[car]
    orientation = get_car_orientation(positions)
    
    # Clear current positions
    for pos in positions:
        board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if orientation == 'H':
            new_positions.append((pos[0], pos[1] + direction))
        else:
            new_positions.append((pos[0] + direction, pos[1]))
    
    # Place car in new positions
    for pos in new_positions:
        board[pos[0]][pos[1]] = car
    
    return [''.join(row) for row in board]

def board_to_string(board):
    return ''.join(board)

def solve_rush_hour():
    initial = [
        'EBBCCC',
        'E....H',
        'F.xAAH',
        'F.G...',
        '..GDDD',
        '......'
    ]
    
    queue = deque([(initial, [])])
    seen = {board_to_string(initial)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if red car (AA) can exit
        red_car = cars['A']
        if max(p[1] for p in red_car) == 4:  # Can move to exit
            return moves
        
        for car, direction in get_valid_moves(current_board, cars):
            new_board = make_move(current_board, cars, car, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format output
solution = solve_rush_hour()
if solution:
    formatted = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted.append(f"{car}{sign}{abs(direction)}")
    print('<<<' + ' '.join(formatted) + '>>>')
else:
    print("No solution found")