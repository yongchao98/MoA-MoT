from collections import deque

# Initial board state
initial_board = [
    list("BBBCCM"),
    list("DDEEEM"),
    list("AAJK.."),
    list("I.JKFF"),
    list("I.GGLx"),
    list("HH..L.")
]

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return cars

def is_path_clear(board, row, start_col, end_col):
    for col in range(start_col, end_col + 1):
        if board[row][col] != '.':
            return False
    return True

def get_blocking_cars(board, red_car_row, red_car_end):
    blocking = []
    for col in range(red_car_end + 1, len(board[0])):
        if board[red_car_row][col] not in ['.', 'x']:
            blocking.append(board[red_car_row][col])
    return blocking

def can_move_vertical(board, col, start_row, end_row):
    for row in range(start_row, end_row + 1):
        if row < 0 or row >= len(board) or board[row][col] != '.':
            return False
    return True

def can_move_horizontal(board, row, start_col, end_col):
    for col in range(start_col, end_col + 1):
        if col < 0 or col >= len(board[0]) or board[row][col] != '.':
            return False
    return True

def get_valid_moves(board, cars):
    moves = []
    red_car_pos = cars['A']
    red_car_row = red_car_pos[0][0]
    red_car_end = max(pos[1] for pos in red_car_pos)
    blocking_cars = get_blocking_cars(board, red_car_row, red_car_end)
    
    for car in blocking_cars + ['A']:
        positions = cars[car]
        is_horizontal = positions[0][0] == positions[-1][0]
        
        if is_horizontal:
            row = positions[0][0]
            left = min(pos[1] for pos in positions)
            right = max(pos[1] for pos in positions)
            
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(pos[0] for pos in positions)
            bottom = max(pos[0] for pos in positions)
            
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, move):
    car, direction = move
    new_board = [row[:] for row in board]
    positions = cars[car]
    is_horizontal = positions[0][0] == positions[-1][0]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    for pos in positions:
        new_row = pos[0]
        new_col = pos[1]
        if is_horizontal:
            new_col += direction
        else:
            new_row += direction
        new_board[new_row][new_col] = car
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        if 'A' in cars:
            red_car_pos = cars['A']
            if max(pos[1] for pos in red_car_pos) == 4:
                return moves
        
        for move in get_valid_moves(current_board, cars):
            new_board = apply_move(current_board, cars, move)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [move]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{car}{sign}1 "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")