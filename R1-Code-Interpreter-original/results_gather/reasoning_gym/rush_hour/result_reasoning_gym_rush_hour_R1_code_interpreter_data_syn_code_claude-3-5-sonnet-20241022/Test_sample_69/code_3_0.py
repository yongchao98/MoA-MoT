from collections import deque

def create_board():
    return [
        list('GBBBI.'),
        list('GCC.I.'),
        list('AAH.IJ'),
        list('.xHDDJ'),
        list('EEE..J'),
        list('..FFF.')
    ]

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[board[i][j]]['positions'].append((i, j))
    
    # Determine orientation and length for each car
    for car, info in cars.items():
        pos = info['positions']
        info['orientation'] = 'H' if pos[0][0] == pos[1][0] else 'V'
        info['length'] = len(pos)
        info['positions'].sort()
    return cars

def can_move(board, car_info, direction):
    pos = car_info['positions']
    if car_info['orientation'] == 'H':
        row = pos[0][0]
        if direction > 0:  # right
            col = max(p[1] for p in pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # left
            col = min(p[1] for p in pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = pos[0][1]
        if direction > 0:  # down
            row = max(p[0] for p in pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # up
            row = min(p[0] for p in pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_info, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for pos in car_info['positions']:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    if car_info['orientation'] == 'H':
        row = car_info['positions'][0][0]
        for pos in car_info['positions']:
            new_col = pos[1] + direction
            new_board[row][new_col] = car
            new_positions.append((row, new_col))
    else:
        col = car_info['positions'][0][1]
        for pos in car_info['positions']:
            new_row = pos[0] + direction
            new_board[new_row][col] = car
            new_positions.append((new_row, col))
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_blocking_red_car(car_info, red_row):
    return (car_info['orientation'] == 'V' and 
            any(pos[0] == red_row and pos[1] > 1 for pos in car_info['positions']))

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_info(current_board)
        
        # Check if solved
        red_car = cars['A']
        if max(pos[1] for pos in red_car['positions']) >= 4:
            return moves
        
        # Get red car row
        red_row = red_car['positions'][0][0]
        
        # First priority: Move blocking cars in red car's path
        for car, info in cars.items():
            if car == 'A':
                continue
                
            if is_blocking_red_car(info, red_row):
                for direction in [-1, 1]:
                    if can_move(current_board, info, direction):
                        new_board = move_car(current_board, car, info, direction)
                        board_str = board_to_string(new_board)
                        
                        if board_str not in seen:
                            seen.add(board_str)
                            move_str = f"{car}{'+' if direction > 0 else '-'}1"
                            new_moves = moves + [move_str]
                            queue.appendleft((new_board, new_moves))
        
        # Second priority: Move red car and other cars
        for car, info in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, info, direction):
                    new_board = move_car(current_board, car, info, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")