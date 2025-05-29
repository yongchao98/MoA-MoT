from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            car = board[i][j]
            if car != '.' and car != 'x':
                if car not in cars:
                    cars[car] = [(i, j)]
                else:
                    cars[car].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction > 0:  # moving right
            col = max(pos[1] for pos in car_positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # moving left
            col = min(pos[1] for pos in car_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # moving down
            row = max(pos[0] for pos in car_positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # moving up
            row = min(pos[0] for pos in car_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        cols = [pos[1] for pos in positions]
        new_cols = [col + direction for col in cols]
        for col in new_cols:
            new_board[row][col] = car
    else:
        col = positions[0][1]
        rows = [pos[0] for pos in positions]
        new_rows = [row + direction for row in rows]
        for row in new_rows:
            new_board[row][col] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (red car AA at exit)
        if 'A' in current_board[2] and current_board[2].rindex('A') == 5:
            return moves
        
        cars = get_car_positions(current_board)
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_tuple = tuple(new_board)
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = [
    'xH.BBB',
    '.H.xKL',
    'AA.JKL',
    'GCCJ.M',
    'G.IDDM',
    'EEIFFF'
]

solution = solve_puzzle(board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")