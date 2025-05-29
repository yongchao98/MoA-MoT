from collections import deque

def create_board():
    return [
        list('BBCC.x'),
        list('DDJEEM'),
        list('..JAAM'),
        list('x.KFFM'),
        list('..KLGG'),
        list('.HHLII')
    ]

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
                car = board[i][j]
                if car not in cars:
                    # Find orientation and length
                    if j + 1 < 6 and board[i][j + 1] == car:  # horizontal
                        length = 2
                        if j + 2 < 6 and board[i][j + 2] == car:
                            length = 3
                        cars[car] = {'pos': (i, j), 'horizontal': True, 'length': length}
                    elif i + 1 < 6 and board[i + 1][j] == car:  # vertical
                        length = 2
                        if i + 2 < 6 and board[i + 2][j] == car:
                            length = 3
                        cars[car] = {'pos': (i, j), 'horizontal': False, 'length': length}
                    else:  # single cell
                        cars[car] = {'pos': (i, j), 'horizontal': True, 'length': 1}
    return cars

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_position(board, i, j):
    return 0 <= i < 6 and 0 <= j < 6 and board[i][j] in ['.', 'x']

def can_move(board, car_info, direction):
    car, info = car_info
    i, j = info['pos']
    length = info['length']
    
    if info['horizontal']:
        if direction > 0:  # right
            new_j = j + length
            return new_j < 6 and board[i][new_j] == '.'
        else:  # left
            return j > 0 and board[i][j-1] == '.' and board[i][j-1] != 'x'
    else:  # vertical
        if direction > 0:  # down
            new_i = i + length
            return new_i < 6 and board[new_i][j] == '.'
        else:  # up
            return i > 0 and board[i-1][j] == '.' and board[i-1][j] != 'x'

def move_car(board, car_info, direction):
    new_board = [row[:] for row in board]
    car, info = car_info
    i, j = info['pos']
    length = info['length']
    
    # Clear current position
    if info['horizontal']:
        for col in range(j, j + length):
            new_board[i][col] = '.'
        # Set new position
        if direction > 0:  # right
            for col in range(j + 1, j + length + 1):
                if new_board[i][col] != 'x':  # Don't overwrite walls
                    new_board[i][col] = car
        else:  # left
            for col in range(j - 1, j + length - 1):
                if new_board[i][col] != 'x':  # Don't overwrite walls
                    new_board[i][col] = car
    else:  # vertical
        for row in range(i, i + length):
            new_board[row][j] = '.'
        # Set new position
        if direction > 0:  # down
            for row in range(i + 1, i + length + 1):
                if new_board[row][j] != 'x':  # Don't overwrite walls
                    new_board[row][j] = car
        else:  # up
            for row in range(i - 1, i + length - 1):
                if new_board[row][j] != 'x':  # Don't overwrite walls
                    new_board[row][j] = car
    
    return new_board

def is_solved(board):
    # Check if red car (AA) can reach the exit
    return board[2][4] == 'A'

def is_valid_board(board):
    # Check if the board state is valid
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
                # Check if car is properly connected
                car = board[i][j]
                if j + 1 < 6 and board[i][j + 1] == car:
                    continue
                if i + 1 < 6 and board[i + 1][j] == car:
                    continue
                if j > 0 and board[i][j - 1] == car:
                    continue
                if i > 0 and board[i - 1][j] == car:
                    continue
    return True

def solve():
    initial_board = create_board()
    visited = {board_to_string(initial_board)}
    queue = deque([(initial_board, [])])
    max_moves = 50  # Prevent infinite loops
    
    while queue:
        current_board, moves = queue.popleft()
        
        if len(moves) > max_moves:
            continue
            
        if is_solved(current_board):
            return moves
        
        cars = get_car_info(current_board)
        for car_info in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, car_info, direction):
                    new_board = move_car(current_board, car_info, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited and is_valid_board(new_board):
                        visited.add(board_str)
                        move = f"{car_info[0]}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move]))
    
    return None

solution = solve()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")