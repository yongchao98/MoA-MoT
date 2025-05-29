from collections import deque

def get_car_info(board):
    # Returns dict of car positions and their orientations (h: horizontal, v: vertical)
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'pos': [(i, j)], 'orient': None}
                else:
                    cars[board[i][j]]['pos'].append((i, j))
                    # Determine orientation
                    if len(cars[board[i][j]]['pos']) == 2:
                        if cars[board[i][j]]['pos'][0][0] == cars[board[i][j]]['pos'][1][0]:
                            cars[board[i][j]]['orient'] = 'h'
                        else:
                            cars[board[i][j]]['orient'] = 'v'
    return cars

def make_move(board, car, direction):
    board = [list(row) for row in board]
    car_info = get_car_info([''.join(row) for row in board])
    
    # Get car positions and orientation
    positions = car_info[car]['pos']
    orientation = car_info[car]['orient']
    
    # Clear current positions
    for i, j in positions:
        board[i][j] = '.'
    
    # Add new positions
    new_positions = []
    for i, j in positions:
        if orientation == 'h':
            new_j = j + direction
            board[i][new_j] = car
        else:
            new_i = i + direction
            board[new_i][j] = car
            
    return [''.join(row) for row in board]

def get_valid_moves(board):
    moves = []
    cars = get_car_info(board)
    
    for car, info in cars.items():
        if info['orient'] == 'h':
            # Try moving left
            leftmost = min(pos[1] for pos in info['pos'])
            if leftmost > 0 and board[info['pos'][0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            # Try moving right
            rightmost = max(pos[1] for pos in info['pos'])
            if rightmost < 5 and board[info['pos'][0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            topmost = min(pos[0] for pos in info['pos'])
            if topmost > 0 and board[topmost-1][info['pos'][0][1]] == '.':
                moves.append((car, -1))
            # Try moving down
            bottommost = max(pos[0] for pos in info['pos'])
            if bottommost < 5 and board[bottommost+1][info['pos'][0][1]] == '.':
                moves.append((car, 1))
    return moves

def solve_puzzle(board):
    start = tuple(board)
    queue = deque([(start, [])])
    seen = {start}
    
    while queue:
        current_state, moves = queue.popleft()
        current_board = list(current_state)
        
        # Check if solved (red car reaches exit)
        cars = get_car_info(current_board)
        if 'A' in cars and max(pos[1] for pos in cars['A']['pos']) == 4:
            return moves
        
        # Get and apply valid moves
        valid_moves = get_valid_moves(current_board)
        for car, direction in valid_moves:
            new_board = make_move(current_board, car, direction)
            new_state = tuple(new_board)
            
            if new_state not in seen:
                seen.add(new_state)
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_state, moves + [move_str]))
    
    return None

# Initial board
initial_board = [
    "G..x..",
    "GBBJ.L",
    "AAIJ.L",
    "CCIDDL",
    ".HEEK.",
    ".HFFK."
]

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")