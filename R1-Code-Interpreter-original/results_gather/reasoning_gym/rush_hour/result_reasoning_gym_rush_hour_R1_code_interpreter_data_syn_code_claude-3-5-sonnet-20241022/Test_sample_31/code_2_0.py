from collections import deque
import copy

def get_car_info(board):
    # Get car positions and orientations
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'pos': [(i, j)], 'orientation': None}
                else:
                    cars[board[i][j]]['pos'].append((i, j))
    
    # Determine orientation for each car
    for car in cars:
        if cars[car]['pos'][0][0] == cars[car]['pos'][1][0]:
            cars[car]['orientation'] = 'H'  # horizontal
        else:
            cars[car]['orientation'] = 'V'  # vertical
    return cars

def is_path_clear(board, row, start_col, end_col):
    for col in range(start_col, end_col):
        if board[row][col] != '.':
            return False
    return True

def solve_rush_hour(board):
    def board_to_string(b):
        return ''.join(''.join(row) for row in b)
    
    def get_moves(current_board, cars, last_car=None):
        moves = []
        for car in cars:
            if car == last_car:  # Avoid moving the same car twice in a row
                continue
            
            positions = cars[car]['pos']
            orientation = cars[car]['orientation']
            
            if orientation == 'H':
                row = positions[0][0]
                left_col = min(p[1] for p in positions)
                right_col = max(p[1] for p in positions)
                
                # Try moving left
                if left_col > 0 and current_board[row][left_col-1] == '.':
                    moves.append((car, -1))
                # Try moving right
                if right_col < 5 and current_board[row][right_col+1] == '.':
                    moves.append((car, 1))
            else:  # vertical
                col = positions[0][1]
                top_row = min(p[0] for p in positions)
                bottom_row = max(p[0] for p in positions)
                
                # Try moving up
                if top_row > 0 and current_board[top_row-1][col] == '.':
                    moves.append((car, -1))
                # Try moving down
                if bottom_row < 5 and current_board[bottom_row+1][col] == '.':
                    moves.append((car, 1))
        return moves

    def apply_move(current_board, cars, move):
        car, direction = move
        new_board = [list(row) for row in current_board]
        positions = cars[car]['pos']
        orientation = cars[car]['orientation']
        
        # Clear current positions
        for pos in positions:
            new_board[pos[0]][pos[1]] = '.'
        
        # Add new positions
        new_positions = []
        for pos in positions:
            if orientation == 'H':
                new_pos = (pos[0], pos[1] + direction)
            else:
                new_pos = (pos[0] + direction, pos[1])
            new_board[new_pos[0]][new_pos[1]] = car
            new_positions.append(new_pos)
        
        return [''.join(row) for row in new_board]

    def get_blocking_cars(board, cars):
        red_car = cars['A']
        red_row = red_car['pos'][0][0]
        red_right = max(p[1] for p in red_car['pos'])
        blocking = []
        
        for col in range(red_right + 1, 5):
            if board[red_row][col] not in '.x':
                blocking.append(board[red_row][col])
        return blocking

    queue = deque([(board, [], None)])
    seen = {board_to_string(board)}
    cars = get_car_info(board)
    
    while queue:
        current_board, path, last_car = queue.popleft()
        cars = get_car_info(current_board)
        
        # Check if solved
        red_car = cars['A']
        if max(p[1] for p in red_car['pos']) == 4:  # Red car reaches the exit
            return path
        
        # Get blocking cars and prioritize their moves
        blocking = get_blocking_cars(current_board, cars)
        moves = get_moves(current_board, cars, last_car)
        
        # Sort moves to prioritize moving blocking cars
        moves.sort(key=lambda x: (x[0] not in blocking))
        
        for move in moves:
            new_board = apply_move(current_board, cars, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, path + [move_str], car))
    
    return None

# Initial board
initial_board = [
    "..xBBL",
    "..ICCL",
    "AAIJ..",
    "HDDJEE",
    "HFFJKx",
    "GG..K."
]

solution = solve_rush_hour(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")