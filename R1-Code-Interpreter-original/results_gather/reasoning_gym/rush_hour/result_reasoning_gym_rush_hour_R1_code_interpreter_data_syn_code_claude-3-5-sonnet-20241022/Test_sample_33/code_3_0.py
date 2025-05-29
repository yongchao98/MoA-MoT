from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n')]

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def move_car(board, car, positions, direction):
    new_board = [row[:] for row in board]
    horizontal = is_horizontal(positions)
    
    # Remove car from old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add car to new positions
    new_positions = []
    for i, j in positions:
        new_i = i + (0 if horizontal else direction)
        new_j = j + (direction if horizontal else 0)
        new_board[new_i][new_j] = car
        new_positions.append((new_i, new_j))
    
    return new_board, new_positions

def can_move(board, positions, direction):
    horizontal = is_horizontal(positions)
    if horizontal:
        row = positions[0][0]
        if direction < 0:  # Move left
            col = min(p[1] for p in positions)
            return col > 0 and board[row][col-1] == '.'
        else:  # Move right
            col = max(p[1] for p in positions)
            return col < 5 and board[row][col+1] == '.'
    else:
        col = positions[0][1]
        if direction < 0:  # Move up
            row = min(p[0] for p in positions)
            return row > 0 and board[row-1][col] == '.'
        else:  # Move down
            row = max(p[0] for p in positions)
            return row < 5 and board[row+1][col] == '.'

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    cars = get_car_positions(board)
    
    # Initial state
    start_state = (board, cars, [])
    visited = {board_to_string(board)}
    queue = deque([start_state])
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if solved (red car at exit)
        red_car_positions = current_cars['A']
        if max(pos[1] for pos in red_car_positions) == 5:
            return moves
        
        # Try all possible moves for each car
        for car, positions in current_cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board, new_positions = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_cars = dict(current_cars)
                        new_cars[car] = new_positions
                        new_moves = moves + [(car, direction)]
                        queue.append((new_board, new_cars, new_moves))
    
    return None

# Initial board
initial_board = """..xBBK
GCCJ.K
GAAJ..
DDIJ..
.HIEE.
.HFFx."""

# Find solution
solution = solve_puzzle(initial_board)

if solution:
    # Format the solution
    move_str = ''
    for car, direction in solution:
        move_str += f'{car}{"+"+str(direction) if direction > 0 else str(direction)} '
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")