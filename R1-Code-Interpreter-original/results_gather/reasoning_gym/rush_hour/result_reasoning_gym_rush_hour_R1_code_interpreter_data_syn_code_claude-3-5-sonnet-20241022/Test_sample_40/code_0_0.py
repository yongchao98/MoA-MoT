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

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    if is_horizontal(positions):
        new_positions = [(p[0], p[1] + direction) for p in positions]
    else:
        new_positions = [(p[0] + direction, p[1]) for p in positions]
    
    # Place car in new positions
    for pos in new_positions:
        new_board[pos[0]][pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        if 'A' in cars:
            a_positions = cars['A']
            if max(p[1] for p in a_positions) == len(current_board[0])-2:
                return moves
        
        # Try all possible moves
        for move in get_valid_moves(current_board, cars):
            new_board = apply_move(current_board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + [f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
initial_board = [
    '.EBBH.',
    '.E.GH.',
    '.AAG..',
    '..FCCC',
    '..F...',
    '.DDD.x'
]

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print('No solution found')