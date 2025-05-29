from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        # Horizontal car
        if is_horizontal(positions):
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if right < 5 and board[row][right+1] == '.':
                moves.append((car, 1))
        # Vertical car
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Try moving down
            if bottom < 5 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, car, direction):
    new_board = [list(row) for row in board]
    positions = cars[car]
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
        
    # Calculate new positions
    if is_horizontal(positions):
        new_positions = [(p[0], p[1] + direction) for p in positions]
    else:
        new_positions = [(p[0] + direction, p[1]) for p in positions]
        
    # Place car in new positions
    for i, j in new_positions:
        new_board[i][j] = car
        
    return [''.join(row) for row in new_board]

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
        
        # Check if solved (red car AA can exit)
        cars = get_car_positions(current_board)
        if any(p[0] == 2 and p[1] == 4 for p in cars['A']):
            return moves
            
        # Try all possible moves
        for car, direction in get_valid_moves(current_board, cars):
            new_board = apply_move(current_board, cars, car, direction)
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