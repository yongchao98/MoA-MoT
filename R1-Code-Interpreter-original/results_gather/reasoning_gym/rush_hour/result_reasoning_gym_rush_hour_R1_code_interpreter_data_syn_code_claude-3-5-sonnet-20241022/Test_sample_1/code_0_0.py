from collections import deque
import copy

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
        # Horizontal car
        if is_horizontal(positions):
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        # Vertical car
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((car, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    
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
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [], [])])  # (board, moves, previous_states)
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves, previous_states = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        a_positions = cars['A']
        if max(p[1] for p in a_positions) == len(current_board[0])-2:
            return moves
        
        # Try all possible moves
        valid_moves = get_valid_moves(current_board, cars)
        for move in valid_moves:
            new_board, new_positions = apply_move(current_board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                car, spaces = move
                move_str = f"{car}{'+' if spaces > 0 else ''}{spaces}"
                queue.append((new_board, moves + [move_str], previous_states + [current_board]))
    
    return None

# Initial board
initial_board = [
    '.F.IBB',
    '.F.I.J',
    '.GAA.J',
    '.GCCCJ',
    '.GHDDD',
    '..H.EE'
]

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print('No solution found')