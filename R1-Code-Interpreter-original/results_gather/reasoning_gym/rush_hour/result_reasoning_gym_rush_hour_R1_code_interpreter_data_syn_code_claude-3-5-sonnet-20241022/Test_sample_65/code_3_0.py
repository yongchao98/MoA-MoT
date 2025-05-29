from collections import deque

def get_vehicle_info(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            char = board[i][j]
            if char not in ['.', 'x']:
                if char not in vehicles:
                    vehicles[char] = {'positions': [(i, j)], 'length': 1}
                else:
                    vehicles[char]['positions'].append((i, j))
                    vehicles[char]['length'] += 1
    
    for v in vehicles.values():
        v['positions'].sort()
        if len(v['positions']) > 1:
            v['horizontal'] = v['positions'][0][0] == v['positions'][1][0]
    return vehicles

def manhattan_distance_to_exit(vehicle_info, board_width):
    if 'A' not in vehicle_info:
        return float('inf')
    positions = vehicle_info['A']['positions']
    rightmost = max(p[1] for p in positions)
    return board_width - 2 - rightmost

def get_blocking_vehicles(board, vehicles):
    if 'A' not in vehicles:
        return set()
    red_car = vehicles['A']
    red_row = red_car['positions'][0][0]
    rightmost = max(p[1] for p in red_car['positions'])
    blocking = set()
    
    # Check all positions from red car to exit
    for col in range(rightmost + 1, len(board[0])):
        if board[red_row][col] not in ['.', 'x']:
            blocking.add(board[red_row][col])
    return blocking

def get_moves(board, vehicles, blocking_vehicles):
    moves = []
    height, width = len(board), len(board[0])
    
    # First try moves for blocking vehicles
    priority_vehicles = list(blocking_vehicles) + ['A']
    other_vehicles = [v for v in vehicles if v not in priority_vehicles]
    all_vehicles = priority_vehicles + other_vehicles
    
    for vehicle in all_vehicles:
        if vehicle not in vehicles:
            continue
        info = vehicles[vehicle]
        positions = info['positions']
        if len(positions) > 1:
            if info['horizontal']:
                leftmost = min(p[1] for p in positions)
                rightmost = max(p[1] for p in positions)
                row = positions[0][0]
                
                if leftmost > 0 and board[row][leftmost-1] == '.':
                    moves.append((vehicle, -1))
                if rightmost < width-1 and board[row][rightmost+1] == '.':
                    moves.append((vehicle, 1))
            else:
                topmost = min(p[0] for p in positions)
                bottommost = max(p[0] for p in positions)
                col = positions[0][1]
                
                if topmost > 0 and board[topmost-1][col] == '.':
                    moves.append((vehicle, -1))
                if bottommost < height-1 and board[bottommost+1][col] == '.':
                    moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    info = vehicles[vehicle]
    new_board = [list(row) for row in board]
    
    for pos in info['positions']:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in info['positions']:
        if info['horizontal']:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    new_vehicles = vehicles.copy()
    new_vehicles[vehicle] = dict(info)
    new_vehicles[vehicle]['positions'] = sorted(new_positions)
    
    return [''.join(row) for row in new_board], new_vehicles

def solve_rush_hour(board):
    vehicles = get_vehicle_info(board)
    start_state = tuple(board)
    visited = {start_state}
    queue = deque([(board, vehicles, [], manhattan_distance_to_exit(vehicles, len(board[0])))])
    
    min_solution = None
    min_length = float('inf')
    
    while queue:
        current_board, current_vehicles, path, distance = queue.popleft()
        
        if len(path) >= min_length:
            continue
            
        # Check if solved
        red_car = current_vehicles['A']
        if red_car['horizontal']:
            rightmost = max(p[1] for p in red_car['positions'])
            if rightmost == len(board[0]) - 2:
                solution = path + [('A', 1)]
                if len(solution) < min_length:
                    min_length = len(solution)
                    min_solution = solution
                continue
        
        blocking = get_blocking_vehicles(current_board, current_vehicles)
        for move in get_moves(current_board, current_vehicles, blocking):
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            new_state = tuple(new_board)
            
            if new_state not in visited:
                visited.add(new_state)
                new_distance = manhattan_distance_to_exit(new_vehicles, len(board[0]))
                queue.append((new_board, new_vehicles, path + [move], new_distance))
    
    return min_solution

board = [
    'E..H..',
    'E..Hxx',
    'EAAIJ.',
    '.FGIJK',
    '.FGBBK',
    'CCDDD.'
]

solution = solve_rush_hour(board)
if solution:
    result = ''
    for vehicle, spaces in solution:
        sign = '+' if spaces > 0 else ''
        result += f'{vehicle}{sign}{spaces} '
    print(f"<<<{result.strip()}>>>")
else:
    print("No solution found")