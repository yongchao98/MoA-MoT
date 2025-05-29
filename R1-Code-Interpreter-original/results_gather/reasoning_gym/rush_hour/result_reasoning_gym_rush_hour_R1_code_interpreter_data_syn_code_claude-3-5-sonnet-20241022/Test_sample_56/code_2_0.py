from collections import deque

def get_blocking_cars(board):
    # Find row of car AA
    aa_row = None
    aa_right = None
    for i, row in enumerate(board):
        if 'A' in row:
            aa_row = i
            aa_right = row.rindex('A')
            break
    
    # Get all cars blocking AA's path to exit
    blocking = set()
    for j in range(aa_right + 1, 6):
        if board[aa_row][j] != '.' and board[aa_row][j] != 'x':
            blocking.add(board[aa_row][j])
    return blocking

def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[board[i][j]]['positions'].append((i, j))
    
    # Determine orientation
    for car in cars:
        pos = cars[car]['positions']
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[-1][0] else 'V'
        cars[car]['positions'].sort()
    
    return cars

def can_move(board, car_info, direction):
    positions = car_info['positions']
    if car_info['orientation'] == 'H':
        row = positions[0][0]
        if direction > 0:  # right
            return (positions[-1][1] + 1 < 6 and 
                   board[row][positions[-1][1] + 1] == '.')
        else:  # left
            return (positions[0][1] - 1 >= 0 and 
                   board[row][positions[0][1] - 1] == '.')
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # down
            return (positions[-1][0] + 1 < 6 and 
                   board[positions[-1][0] + 1][col] == '.')
        else:  # up
            return (positions[0][0] - 1 >= 0 and 
                   board[positions[0][0] - 1][col] == '.')

def move_car(board, car, car_info, direction):
    new_board = [list(row) for row in board]
    positions = car_info['positions']
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if car_info['orientation'] == 'H':
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board]

def board_to_string(board):
    return '\n'.join(board)

def solve_puzzle():
    initial_board = [
        "..HBBB",
        "..HICC",
        ".AAIJK",
        ".GDDJK",
        ".GEEJL",
        "FF...L"
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    
    while queue:
        current_board, moves = queue.popleft()
        board_str = board_to_string(current_board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        # Check if solved
        for row in current_board:
            if 'A' in row and row.rindex('A') == 5:
                return moves
        
        cars = get_car_info(current_board)
        blocking = get_blocking_cars(current_board)
        
        # Prioritize moves for blocking cars and the red car
        priority_cars = ['A'] + list(blocking)
        other_cars = [c for c in cars if c not in priority_cars]
        all_cars = priority_cars + other_cars
        
        for car in all_cars:
            car_info = cars[car]
            
            # Try both directions
            for direction in [-1, 1]:
                if can_move(current_board, car_info, direction):
                    new_board = move_car(current_board, car, car_info, direction)
                    move_str = f"{car}{'+' if direction > 0 else '-'}1"
                    
                    # Prioritize moves that clear path for red car
                    new_blocking = get_blocking_cars(new_board)
                    if car == 'A' or len(new_blocking) < len(blocking):
                        queue.appendleft((new_board, moves + [move_str]))
                    else:
                        queue.append((new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")