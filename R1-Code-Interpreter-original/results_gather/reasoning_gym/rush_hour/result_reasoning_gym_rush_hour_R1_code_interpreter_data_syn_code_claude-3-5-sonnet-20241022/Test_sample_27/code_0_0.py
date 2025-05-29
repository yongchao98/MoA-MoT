from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            car = board[i][j]
            if car != '.' and car != 'x':
                if car not in cars:
                    cars[car] = []
                cars[car].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_possible_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Try moving down
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, car, direction):
    new_board = [row[:] for row in board]
    new_cars = copy.deepcopy(cars)
    positions = new_cars[car]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Update positions
    new_positions = []
    for i, j in positions:
        if horizontal:
            new_positions.append((i, j + direction))
        else:
            new_positions.append((i + direction, j))
    
    # Set new positions
    for i, j in new_positions:
        new_board[i][j] = car
    new_cars[car] = new_positions
    
    return new_board, new_cars

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    cars = get_car_positions(board)
    
    queue = deque([(board, cars, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if red car (AA) has reached the exit
        red_car_positions = current_cars['A']
        if max(p[1] for p in red_car_positions) == len(board[0])-1:
            return moves
        
        possible_moves = get_possible_moves(current_board, current_cars)
        for car, direction in possible_moves:
            new_board, new_cars = apply_move(current_board, current_cars, car, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_cars, new_moves))
    
    return None

# Initial board
board = """FBB.J.
F.I.JK
AAI..K
CCI..L
GHDDDL
GHEEx."""

solution = solve_puzzle(board)
if solution:
    # Format the solution
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("No solution found")