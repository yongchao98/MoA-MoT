from collections import deque
import copy

# Initial board state
initial_board = [
    list("BBBCCM"),
    list("DDEEEM"),
    list("AAJK.."),
    list("I.JKFF"),
    list("I.GGLx"),
    list("HH..L.")
]

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
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
    positions = cars[car]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        new_pos = (pos[0], pos[1] + direction) if horizontal else (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_cars = get_car_positions(initial_board)
    queue = deque([(initial_board, [], initial_cars)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, cars = queue.popleft()
        
        # Check if red car (AA) can exit
        if 'A' in cars and cars['A'][0][0] == 2 and max(p[1] for p in cars['A']) == 4:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars)
        for car, direction in valid_moves:
            new_board = apply_move(current_board, cars, car, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_cars = get_car_positions(new_board)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves, new_cars))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{car}{sign}1 "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")