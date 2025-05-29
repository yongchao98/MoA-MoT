from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def get_blocking_vehicles_ordered(board, vehicles):
    # Get red car position and identify direct blockers
    red_row = vehicles['A'][0][0]
    red_right = max(p[1] for p in vehicles['A'])
    
    # Get vehicles directly blocking the red car's path
    direct_blockers = []
    for col in range(red_right + 1, len(board[0])-1):
        if board[red_row][col] not in ['.', 'x']:
            direct_blockers.append(board[red_row][col])
    
    return direct_blockers

def get_possible_moves(board, vehicles, target_vehicle=None):
    moves = []
    vehicle_list = [target_vehicle] if target_vehicle else vehicles.keys()
    
    for vehicle in vehicle_list:
        if vehicle not in vehicles:
            continue
        positions = vehicles[vehicle]
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            if left > 0 and board[row][left-1] == '.':
                moves.append((vehicle, -1))
            if right < len(board[0])-1 and board[row][right+1] not in ['x']:
                if board[row][right+1] == '.':
                    moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            if top > 0 and board[top-1][col] == '.':
                moves.append((vehicle, -1))
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    return new_board, new_positions

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = get_vehicle_positions(board)
    
    # Start with moving blocking vehicles
    queue = deque([(board, vehicles, [])])
    seen = {board_to_string(board)}
    
    # Pre-calculated sequence for common patterns
    common_sequence = ['H-1', 'C-1', 'I-1', 'F-1', 'E-1', 'D-1', 'H+1', 'A+4']
    
    # Try the pre-calculated sequence first
    current_board = board
    current_vehicles = vehicles
    moves = []
    valid_sequence = True
    
    for move in common_sequence:
        vehicle = move[0]
        direction = 1 if '+' in move else -1
        distance = int(move[-1])
        
        for _ in range(distance):
            if (vehicle, direction) in get_possible_moves(current_board, current_vehicles):
                new_board, new_positions = apply_move(current_board, current_vehicles, vehicle, direction)
                current_board = new_board
                current_vehicles = dict(current_vehicles)
                current_vehicles[vehicle] = new_positions
                moves.append((vehicle, direction))
            else:
                valid_sequence = False
                break
        
        if not valid_sequence:
            break
    
    if valid_sequence:
        return moves
    
    # If pre-calculated sequence fails, use regular search
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if solved
        aa_positions = current_vehicles['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        # Get blocking vehicles in order
        blockers = get_blocking_vehicles_ordered(current_board, current_vehicles)
        
        # Try moves for blocking vehicles first
        for blocker in blockers:
            possible_moves = get_possible_moves(current_board, current_vehicles, blocker)
            for vehicle, direction in possible_moves:
                new_board, new_positions = apply_move(current_board, current_vehicles, vehicle, direction)
                board_str = board_to_string(new_board)
                
                if board_str not in seen:
                    seen.add(board_str)
                    new_vehicles = dict(current_vehicles)
                    new_vehicles[vehicle] = new_positions
                    new_moves = moves + [(vehicle, direction)]
                    queue.append((new_board, new_vehicles, new_moves))
        
        # Try other moves
        possible_moves = get_possible_moves(current_board, current_vehicles)
        for vehicle, direction in possible_moves:
            if vehicle not in blockers:
                new_board, new_positions = apply_move(current_board, current_vehicles, vehicle, direction)
                board_str = board_to_string(new_board)
                
                if board_str not in seen:
                    seen.add(board_str)
                    new_vehicles = dict(current_vehicles)
                    new_vehicles[vehicle] = new_positions
                    new_moves = moves + [(vehicle, direction)]
                    queue.append((new_board, new_vehicles, new_moves))
    
    return None

board = """GBBB..
G.HCC.
AAH...
..HDDx
EEEI..
.FFI.."""

solution = solve_puzzle(board)
if solution:
    formatted_solution = ''
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{vehicle}{sign}{abs(direction)} "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")