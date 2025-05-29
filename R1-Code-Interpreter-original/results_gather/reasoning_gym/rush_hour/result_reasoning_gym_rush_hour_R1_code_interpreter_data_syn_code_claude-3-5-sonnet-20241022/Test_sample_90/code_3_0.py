from collections import deque

def create_board():
    return [
        list("BBBKCC"),
        list("DDJK.L"),
        list("IJAAL"),
        list("IEE.xM"),
        list("FF...M"),
        list("GGHH.x")
    ]

def find_cars(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            car = board[i][j]
            if car not in '.x':
                if car not in cars:
                    cars[car] = {'pos': [(i, j)]}
                else:
                    cars[car]['pos'].append((i, j))
    
    # Determine orientation
    for car in cars:
        pos = cars[car]['pos']
        cars[car]['horizontal'] = pos[0][0] == pos[1][0]
    return cars

def move_car(board, car_id, car_data, delta):
    new_board = [row[:] for row in board]
    new_pos = []
    
    # Remove car from old position
    for i, j in car_data['pos']:
        new_board[i][j] = '.'
    
    # Add car to new position
    for i, j in car_data['pos']:
        if car_data['horizontal']:
            new_j = j + delta
            new_board[i][new_j] = car_id
            new_pos.append((i, new_j))
        else:
            new_i = i + delta
            new_board[new_i][j] = car_id
            new_pos.append((new_i, j))
    
    return new_board, new_pos

def is_valid_move(board, car_data, delta):
    if car_data['horizontal']:
        row = car_data['pos'][0][0]
        if delta > 0:  # Moving right
            col = max(p[1] for p in car_data['pos']) + delta
            return col < 6 and board[row][col] == '.'
        else:  # Moving left
            col = min(p[1] for p in car_data['pos']) + delta
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_data['pos'][0][1]
        if delta > 0:  # Moving down
            row = max(p[0] for p in car_data['pos']) + delta
            return row < 6 and board[row][col] == '.'
        else:  # Moving up
            row = min(p[0] for p in car_data['pos']) + delta
            return row >= 0 and board[row][col] == '.'

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve():
    initial_board = create_board()
    cars = find_cars(initial_board)
    
    queue = deque([(initial_board, cars, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if solved (red car at exit)
        red_car = current_cars['A']
        if max(p[1] for p in red_car['pos']) == 4:
            return moves
        
        # Try all possible moves
        for car_id, car_data in current_cars.items():
            for delta in [-1, 1]:
                if is_valid_move(current_board, car_data, delta):
                    new_board, new_pos = move_car(current_board, car_id, car_data, delta)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_cars = dict(current_cars)
                        new_cars[car_id] = {'pos': new_pos, 'horizontal': car_data['horizontal']}
                        move = f"{car_id}{'+' if delta > 0 else '-'}1"
                        queue.append((new_board, new_cars, moves + [move]))
    
    return None

solution = solve()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print('<<<No solution found>>>')