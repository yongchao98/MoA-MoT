from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

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
        # Horizontal car
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
        # Vertical car
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
    new_positions = []
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    if is_horizontal(positions):
        for pos in positions:
            new_pos = (pos[0], pos[1] + direction)
            new_positions.append(new_pos)
            new_board[new_pos[0]][new_pos[1]] = car
    else:
        for pos in positions:
            new_pos = (pos[0] + direction, pos[1])
            new_positions.append(new_pos)
            new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle(initial_board):
    board = initial_board.split('\n')
    visited = set()
    queue = deque([(board, [], {})])
    
    while queue:
        current_board, moves, _ = queue.popleft()
        board_str = '\n'.join(current_board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        cars = get_car_positions(current_board)
        
        # Check if solved (AA reaches right edge)
        aa_positions = cars['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        # Try all possible moves
        valid_moves = get_valid_moves(current_board, cars)
        for car, direction in valid_moves:
            new_board, new_positions = apply_move(current_board, cars, car, direction)
            new_moves = moves + [(car, direction)]
            queue.append((new_board, new_moves, {car: new_positions}))
    
    return None

# Initial board
board = """BBxCCK
G..I.K
GAAI.L
..HDDL
..H.J.
EEFFJx"""

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