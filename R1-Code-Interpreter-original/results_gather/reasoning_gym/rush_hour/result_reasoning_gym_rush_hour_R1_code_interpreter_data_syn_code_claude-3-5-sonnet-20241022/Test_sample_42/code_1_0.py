from collections import deque

def get_blocking_cars(board):
    # Find cars blocking the path of the red car (AA)
    blocking = set()
    red_car_row = None
    for i in range(6):
        if 'A' in board[i]:
            red_car_row = i
            break
    
    if red_car_row is not None:
        # Find the rightmost position of AA
        right_pos = max(j for j, c in enumerate(board[red_car_row]) if c == 'A')
        # Check all positions to the right of AA
        for j in range(right_pos + 1, 6):
            if board[red_car_row][j] not in '.x':
                blocking.add(board[red_car_row][j])
    return blocking

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

def get_valid_moves(board, cars, blocking_cars):
    moves = []
    priority_moves = []  # Moves for blocking cars
    
    for car, positions in cars.items():
        if len(positions) < 2:
            continue
            
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                move = (car, -1)
                if car in blocking_cars:
                    priority_moves.append(move)
                else:
                    moves.append(move)
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < 5 and board[positions[0][0]][rightmost+1] == '.':
                move = (car, 1)
                if car in blocking_cars:
                    priority_moves.append(move)
                else:
                    moves.append(move)
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                move = (car, -1)
                if car in blocking_cars:
                    priority_moves.append(move)
                else:
                    moves.append(move)
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < 5 and board[bottommost+1][positions[0][1]] == '.':
                move = (car, 1)
                if car in blocking_cars:
                    priority_moves.append(move)
                else:
                    moves.append(move)
    
    return priority_moves + moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        blocking_cars = get_blocking_cars(current_board)
        
        # Check if solved
        if 'A' in cars and max(p[1] for p in cars['A']) == 4:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars, blocking_cars)
        for move in valid_moves:
            new_board = apply_move(current_board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, moves + [move_str]))
    
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