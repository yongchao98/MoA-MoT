from heapq import heappush, heappop

def create_board():
    return [
        ['.', '.', 'x', 'B', 'B', 'B'],
        ['G', 'H', 'I', '.', 'C', 'C'],
        ['G', 'H', 'I', 'A', 'A', 'K'],
        ['D', 'D', 'D', 'J', '.', 'K'],
        ['.', '.', '.', 'J', 'E', 'E'],
        ['.', 'F', 'F', 'F', '.', 'x']
    ]

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'pos': [], 'orientation': None}
                cars[board[i][j]]['pos'].append((i, j))
    
    for car in cars:
        pos = cars[car]['pos']
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[1][0] else 'V'
        cars[car]['pos'].sort()
    return cars

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board, car_info):
    # Distance of red car from exit
    red_car = car_info['A']
    red_pos = red_car['pos'][-1]
    exit_pos = (2, 5)
    
    # Count blocking vehicles
    blocking = 0
    for car, info in car_info.items():
        if car != 'A':
            for pos in info['pos']:
                if pos[0] == 2 and pos[1] > red_pos[1]:
                    blocking += 1
                    break
    
    return manhattan_distance(red_pos, exit_pos) + blocking * 2

def make_move(board, car, car_info, direction):
    new_board = [row[:] for row in board]
    positions = car_info[car]['pos']
    orientation = car_info[car]['orientation']
    
    if orientation == 'H':
        row = positions[0][0]
        start_col = positions[0][1]
        length = len(positions)
        new_col = start_col + direction
        
        if new_col < 0 or new_col + length > 6:
            return None
            
        # Check if move is valid
        if direction > 0 and (new_col + length > 6 or new_board[row][new_col + length - 1] not in '.'):
            return None
        if direction < 0 and (new_col < 0 or new_board[row][new_col] not in '.'):
            return None
            
        # Make move
        for j in range(len(positions)):
            new_board[row][positions[j][1]] = '.'
        for j in range(length):
            new_board[row][new_col + j] = car
            
    else:  # Vertical
        col = positions[0][1]
        start_row = positions[0][0]
        length = len(positions)
        new_row = start_row + direction
        
        if new_row < 0 or new_row + length > 6:
            return None
            
        # Check if move is valid
        if direction > 0 and (new_row + length > 6 or new_board[new_row + length - 1][col] not in '.'):
            return None
        if direction < 0 and (new_row < 0 or new_board[new_row][col] not in '.'):
            return None
            
        # Make move
        for i in range(len(positions)):
            new_board[positions[i][0]][col] = '.'
        for i in range(length):
            new_board[new_row + i][col] = car
            
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    initial_state = board_to_string(initial_board)
    visited = {initial_state}
    pq = [(0, 0, initial_board, [])]  # (priority, moves_made, board, moves)
    
    while pq:
        _, moves_made, current_board, moves = heappop(pq)
        car_info = get_car_info(current_board)
        
        # Check if solved
        red_car = car_info['A']
        red_pos = red_car['pos'][-1]
        if red_pos[1] < 5:  # If not at exit
            can_exit = True
            for j in range(red_pos[1] + 1, 6):
                if current_board[2][j] not in '.x':
                    can_exit = False
                    break
            if can_exit:
                return moves + [f"A+{5 - red_pos[1]}"]
        
        # Try moves for each car
        for car in car_info:
            for direction in [-1, 1]:
                new_board = make_move(current_board, car, car_info, direction)
                if new_board:
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}1"]
                        h = heuristic(new_board, get_car_info(new_board))
                        priority = moves_made + 1 + h
                        heappush(pq, (priority, moves_made + 1, new_board, new_moves))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')