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
    
    # Determine orientation for each vehicle
    for v in vehicles.values():
        pos = v['positions']
        v['horizontal'] = pos[0][0] == pos[1][0]  # Same row = horizontal
        v['positions'].sort()  # Sort positions for consistency
    return vehicles

def make_board_state(board):
    return tuple(tuple(row) for row in board)

def get_moves(board, vehicles):
    moves = []
    for vehicle, info in vehicles.items():
        if info['horizontal']:
            # Get leftmost and rightmost positions
            left_pos = min(p[1] for p in info['positions'])
            right_pos = max(p[1] for p in info['positions'])
            row = info['positions'][0][0]
            
            # Check left move
            if left_pos > 0 and board[row][left_pos - 1] == '.':
                moves.append((vehicle, -1))
            
            # Check right move
            if right_pos < len(board[0]) - 1 and board[row][right_pos + 1] == '.':
                moves.append((vehicle, 1))
        else:
            # Get topmost and bottommost positions
            top_pos = min(p[0] for p in info['positions'])
            bottom_pos = max(p[0] for p in info['positions'])
            col = info['positions'][0][1]
            
            # Check up move
            if top_pos > 0 and board[top_pos - 1][col] == '.':
                moves.append((vehicle, -1))
            
            # Check down move
            if bottom_pos < len(board) - 1 and board[bottom_pos + 1][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    info = vehicles[vehicle]
    new_board = [list(row) for row in board]
    
    # Clear current positions
    for pos in info['positions']:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in info['positions']:
        if info['horizontal']:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    # Update vehicle info
    new_vehicles = vehicles.copy()
    new_vehicles[vehicle] = dict(info)
    new_vehicles[vehicle]['positions'] = sorted(new_positions)
    
    return [''.join(row) for row in new_board], new_vehicles

def solve_rush_hour(board):
    vehicles = get_vehicle_info(board)
    start_state = make_board_state(board)
    visited = {start_state}
    queue = deque([(board, vehicles, [])])
    
    while queue:
        current_board, current_vehicles, path = queue.popleft()
        
        # Check if solved (red car 'A' can reach exit)
        red_car = current_vehicles['A']
        if red_car['horizontal']:
            rightmost = max(p[1] for p in red_car['positions'])
            if rightmost == len(board[0]) - 2:  # One move away from exit
                return path + [('A', 1)]
        
        # Try all possible moves
        for move in get_moves(current_board, current_vehicles):
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            new_state = make_board_state(new_board)
            
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_board, new_vehicles, path + [move]))
    
    return None

# Initial board
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
    # Format the solution
    result = ''
    for vehicle, spaces in solution:
        sign = '+' if spaces > 0 else ''
        result += f'{vehicle}{sign}{spaces} '
    print(f"<<<{result.strip()}>>>")
else:
    print("No solution found")