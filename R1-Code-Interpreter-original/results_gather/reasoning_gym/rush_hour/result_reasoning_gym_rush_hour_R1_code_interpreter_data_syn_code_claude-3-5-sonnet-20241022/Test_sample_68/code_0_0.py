from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
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
        # Horizontal movement
        if is_horizontal(positions):
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            
            # Try moving right
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        
        # Vertical movement
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
    new_board = [list(row) for row in board]
    positions = cars[car]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    if is_horizontal(positions):
        for pos in positions:
            new_positions.append((pos[0], pos[1] + direction))
    else:
        for pos in positions:
            new_positions.append((pos[0] + direction, pos[1]))
    
    # Place car in new positions
    for pos in new_positions:
        new_board[pos[0]][pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "E.FBBJ",
        "E.F.IJ",
        "EAAGIK",
        "...G.K",
        "CC.H.K",
        "DDDH.."
    ]
    
    target_x = 5  # Exit position
    target_y = 2  # Row where AA is
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved
        if 'A' in cars:
            aa_positions = cars['A']
            if max(p[1] for p in aa_positions) == len(current_board[0])-2:
                return moves + [('A', 1)]
        
        valid_moves = get_valid_moves(current_board, cars)
        
        for car, direction in valid_moves:
            new_board = apply_move(current_board, cars, car, direction)
            board_tuple = tuple(new_board)
            
            if board_tuple not in visited:
                visited.add(board_tuple)
                queue.append((new_board, moves + [(car, direction)]))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{car}{sign}{abs(direction)} "
    print(formatted_solution.strip())
else:
    print("No solution found")