from heapq import heappush, heappop

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

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board):
    cars = {}
    red_car_pos = None
    blocking_cars = 0
    
    # Find positions of all cars
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                red_car_pos = (i, j)
            elif board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    
    # Count blocking cars
    red_car_row = red_car_pos[0]
    red_car_rightmost = max(pos[1] for pos in cars['A'])
    
    for j in range(red_car_rightmost + 1, len(board[0])):
        if board[red_car_row][j] not in '.x':
            blocking_cars += 1
    
    return blocking_cars * 2

def solve_puzzle():
    initial_board = [
        'GBB.KL',
        'GHI.KL',
        'GHIAAM',
        'CCCJ.M',
        '..xJDD',
        'EEFF..'
    ]
    
    visited = set()
    pq = [(0, 0, initial_board, [])]  # (priority, moves_count, board, moves)
    visited.add(tuple(initial_board))
    
    while pq:
        _, moves_count, current_board, moves = heappop(pq)
        
        # Check if solved
        cars = get_car_positions(current_board)
        if max(pos[1] for pos in cars['A']) == len(current_board[0])-1:
            return moves
        
        # Try moving each car
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_tuple = tuple(new_board)
                    
                    if board_tuple not in visited:
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        visited.add(board_tuple)
                        h_score = heuristic(new_board)
                        priority = moves_count + 1 + h_score
                        heappush(pq, (priority, moves_count + 1, new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")