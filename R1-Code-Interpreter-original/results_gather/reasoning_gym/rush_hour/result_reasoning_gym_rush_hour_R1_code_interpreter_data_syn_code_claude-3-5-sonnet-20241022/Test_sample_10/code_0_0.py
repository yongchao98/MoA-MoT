from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
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
    if is_horizontal(positions):
        row = positions[0][0]
        for col in [p[1] for p in positions]:
            new_board[row][col] = '.'
        for col in [p[1] + direction for p in positions]:
            new_board[row][col] = car
    else:
        col = positions[0][1]
        for row in [p[0] for p in positions]:
            new_board[row][col] = '.'
        for row in [p[0] + direction for p in positions]:
            new_board[row][col] = car
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        '.xBBCC',
        '..x.K.',
        'G.AAK.',
        'G.IJDD',
        'H.IJ..',
        'HEEFFF'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {'\n'.join(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (AA is at the exit)
        aa_positions = cars['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        # Try all possible moves
        for car, direction in get_valid_moves(current_board, cars):
            new_board = apply_move(current_board, cars, car, direction)
            board_key = '\n'.join(new_board)
            
            if board_key not in seen:
                seen.add(board_key)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f'{car}{sign}1 '
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")