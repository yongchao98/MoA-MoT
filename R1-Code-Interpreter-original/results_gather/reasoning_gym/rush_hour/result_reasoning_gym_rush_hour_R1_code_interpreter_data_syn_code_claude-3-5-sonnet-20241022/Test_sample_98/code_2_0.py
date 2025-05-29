from collections import deque

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
        if direction > 0:
            col = max(pos[1] for pos in car_positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:
            col = min(pos[1] for pos in car_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_positions[0][1]
        if direction > 0:
            row = max(pos[0] for pos in car_positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:
            row = min(pos[0] for pos in car_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    for row, col in positions:
        new_board[row][col] = '.'
    
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

def blocking_cars(board, target_row=2):
    # Find cars blocking the path to the exit for car A
    blocking = set()
    a_pos = board[target_row].index('A')
    for j in range(a_pos + 2, len(board[target_row])):
        if board[target_row][j] != '.' and board[target_row][j] != 'x':
            blocking.add(board[target_row][j])
    return blocking

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        if 'A' in current_board[2] and current_board[2].rindex('A') == 5:
            return moves
        
        # Get blocking cars and prioritize moving them
        blockers = blocking_cars(current_board)
        cars = get_car_positions(current_board)
        
        # First try to move blocking cars
        for car in blockers:
            positions = cars[car]
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_tuple = tuple(new_board)
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        queue.appendleft((new_board, moves + [move_str]))
        
        # Then try to move car A if possible
        if 'A' in cars:
            positions = cars['A']
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_tuple = tuple(new_board)
                    if board_tuple not in seen:
                        seen.add(board_tuple)
                        move_str = f"A{'+' if direction > 0 else '-'}1"
                        queue.appendleft((new_board, moves + [move_str]))
        
        # Finally try other cars
        for car, positions in cars.items():
            if car not in blockers and car != 'A':
                for direction in [-1, 1]:
                    if can_move(current_board, positions, direction):
                        new_board = move_car(current_board, car, positions, direction)
                        board_tuple = tuple(new_board)
                        if board_tuple not in seen:
                            seen.add(board_tuple)
                            move_str = f"{car}{'+' if direction > 0 else '-'}1"
                            queue.append((new_board, moves + [move_str]))
    
    return None

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