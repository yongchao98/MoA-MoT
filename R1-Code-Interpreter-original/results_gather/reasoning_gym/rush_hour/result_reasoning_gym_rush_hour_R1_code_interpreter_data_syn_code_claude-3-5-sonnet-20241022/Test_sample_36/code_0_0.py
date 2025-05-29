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
        if len(positions) < 2:
            continue
        
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
    new_board = [list(row) for row in board]
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
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle(initial_board):
    board = [list(row) for row in initial_board]
    visited = set([''.join(row) for row in board])
    queue = deque([(board, [], {})])
    
    while queue:
        current_board, moves, prev_cars = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car AA reaches right edge)
        aa_positions = cars['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars)
        for car, direction in valid_moves:
            new_board, new_positions = apply_move(current_board, cars, car, direction)
            board_str = ''.join(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves, cars))
    
    return None

# Initial board
initial_board = [
    "BBH.J.",
    "..HxJK",
    "AAH..K",
    "GCCC.K",
    "GDDI..",
    "EE.IFF"
]

# Solve the puzzle
solution = solve_puzzle(initial_board)

# Format the solution
if solution:
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")