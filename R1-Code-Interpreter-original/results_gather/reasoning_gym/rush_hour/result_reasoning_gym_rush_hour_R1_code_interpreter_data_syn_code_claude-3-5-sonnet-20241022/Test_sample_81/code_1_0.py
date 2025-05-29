from collections import deque

def get_board_state():
    return [
        ['.', '.', '.', 'E', '.', '.'],
        ['.', '.', 'C', 'E', '.', '.'],
        ['A', 'A', 'C', 'E', '.', '.'],
        ['.', '.', '.', '.', '.', '.'],
        ['.', '.', 'D', '.', '.', '.'],
        ['.', '.', 'D', 'B', 'B', 'x']
    ]

def is_valid_pos(i, j):
    return 0 <= i < 6 and 0 <= j < 6

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                car = board[i][j]
                if car not in cars:
                    cars[car] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[car]['positions'].append((i, j))
    
    # Determine orientation
    for car in cars:
        pos = cars[car]['positions']
        cars[car]['orientation'] = 'vertical' if len(pos) == 1 or pos[0][1] == pos[-1][1] else 'horizontal'
    
    return cars

def make_move(board, car, direction):
    new_board = [row[:] for row in board]
    cars = get_car_info(board)
    
    if car not in cars:
        return None
    
    car_data = cars[car]
    positions = car_data['positions']
    orientation = car_data['orientation']
    
    # Validate move direction against orientation
    if (orientation == 'horizontal' and direction in ['up', 'down']) or \
       (orientation == 'vertical' and direction in ['left', 'right']):
        return None
    
    # Calculate new positions
    new_positions = []
    for i, j in positions:
        new_i, new_j = i, j
        if direction == 'up': new_i -= 1
        elif direction == 'down': new_i += 1
        elif direction == 'left': new_j -= 1
        elif direction == 'right': new_j += 1
        
        if not is_valid_pos(new_i, new_j) or \
           (new_board[new_i][new_j] not in ['.'] and new_board[new_i][new_j] != car):
            return None
        new_positions.append((new_i, new_j))
    
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in new_positions:
        if new_board[i][j] == 'x':  # Don't overwrite walls
            return None
        new_board[i][j] = car
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = get_board_state()
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_info(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        a_positions = cars['A']['positions']
        if any(j == 4 for i, j in a_positions):  # Check if A can reach exit
            return moves
        
        # Try all possible moves for each car
        for car in cars:
            for direction in ['up', 'down', 'left', 'right']:
                new_board = make_move(current_board, car, direction)
                if new_board:
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move_notation = f"{car}{'+' if direction in ['right', 'down'] else '-'}1"
                        queue.append((new_board, moves + [move_notation]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("<<<No solution found>>>")