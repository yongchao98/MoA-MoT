from collections import deque

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

def get_car_orientation(positions):
    if positions[0][0] == positions[-1][0]:  # same row
        return 'horizontal'
    return 'vertical'

def is_valid_move(board, car_pos, direction, steps):
    orientation = get_car_orientation(car_pos)
    
    if orientation == 'horizontal':
        row = car_pos[0][0]
        if direction > 0:  # moving right
            end_pos = car_pos[-1][1] + steps
            if end_pos >= 6:
                return False
            for i in range(1, steps + 1):
                if car_pos[-1][1] + i >= 6 or board[row][car_pos[-1][1] + i] != '.':
                    return False
            return True
        else:  # moving left
            start_pos = car_pos[0][1] + direction
            if start_pos < 0:
                return False
            for i in range(1, abs(steps) + 1):
                if car_pos[0][1] - i < 0 or board[row][car_pos[0][1] - i] != '.':
                    return False
            return True
    else:  # vertical
        col = car_pos[0][1]
        if direction > 0:  # moving down
            end_pos = car_pos[-1][0] + steps
            if end_pos >= 6:
                return False
            for i in range(1, steps + 1):
                if car_pos[-1][0] + i >= 6 or board[car_pos[-1][0] + i][col] != '.':
                    return False
            return True
        else:  # moving up
            start_pos = car_pos[0][0] + direction
            if start_pos < 0:
                return False
            for i in range(1, abs(steps) + 1):
                if car_pos[0][0] - i < 0 or board[car_pos[0][0] - i][col] != '.':
                    return False
            return True

def make_move(board, car, car_pos, direction):
    new_board = [list(row) for row in board]
    orientation = get_car_orientation(car_pos)
    
    # Clear current position
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new position
    if orientation == 'horizontal':
        for pos in car_pos:
            new_board[pos[0]][pos[1] + direction] = car
    else:
        for pos in car_pos:
            new_board[pos[0] + direction][pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (AA reaches right edge)
        if 'A' in cars and cars['A'][-1][1] == 4:  # Position before exit
            return moves
        
        # Try all possible moves for each car
        for car, positions in cars.items():
            # Try moving 1-3 spaces (reduced from 5 for efficiency)
            for steps in range(1, 4):
                for direction in [-1, 1]:
                    if is_valid_move(current_board, positions, direction, steps):
                        new_board = make_move(current_board, car, positions, direction * steps)
                        new_board_tuple = tuple(new_board)
                        
                        if new_board_tuple not in seen:
                            seen.add(new_board_tuple)
                            move_notation = f"{car}{'+' if direction * steps > 0 else ''}{direction * steps}"
                            queue.append((new_board, moves + [move_notation]))
    
    return None

# Initial board
board = [
    "GBBB..",
    "G..IJ.",
    "HAAIJ.",
    "H.....",
    "HxCCDD",
    "EEFF.."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")