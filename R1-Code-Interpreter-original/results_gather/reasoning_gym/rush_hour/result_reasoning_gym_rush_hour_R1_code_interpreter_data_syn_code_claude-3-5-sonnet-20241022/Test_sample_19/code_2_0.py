from collections import deque, defaultdict

def get_vehicle_positions(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                vehicles[board[i][j]].append((i, j))
    return dict(vehicles)

def get_blocking_vehicles(board, vehicles):
    # Get vehicles directly blocking AA's path
    if 'A' not in vehicles:
        return set()
    aa_row = vehicles['A'][0][0]
    aa_right = vehicles['A'][-1][1]
    blocking = set()
    for col in range(aa_right + 1, len(board[0])):
        for vehicle, positions in vehicles.items():
            if any(pos[0] == aa_row and pos[1] == col for pos in positions):
                blocking.add(vehicle)
    return blocking

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, vehicle_pos, direction):
    if is_horizontal(vehicle_pos):
        if direction > 0:
            return (vehicle_pos[-1][1] + 1 < len(board[0]) and 
                   board[vehicle_pos[0][0]][vehicle_pos[-1][1] + 1] == '.')
        else:
            return (vehicle_pos[0][1] - 1 >= 0 and 
                   board[vehicle_pos[0][0]][vehicle_pos[0][1] - 1] == '.')
    else:
        if direction > 0:
            return (vehicle_pos[-1][0] + 1 < len(board) and 
                   board[vehicle_pos[-1][0] + 1][vehicle_pos[0][1]] == '.')
        else:
            return (vehicle_pos[0][0] - 1 >= 0 and 
                   board[vehicle_pos[0][0] - 1][vehicle_pos[0][1]] == '.')

def move_vehicle(board, vehicle, positions, direction):
    new_board = [list(row) for row in board]
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        for pos in positions:
            new_board[pos[0]][pos[1] + direction] = vehicle
    else:
        for pos in positions:
            new_board[pos[0] + direction][pos[1]] = vehicle
    
    return new_board

def solve_puzzle():
    initial_board = [
        list('.FBBJK'),
        list('.F.IJK'),
        list('xAAI.L'),
        list('CCHI.L'),
        list('EGH...'),
        list('EGDDD.')
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if 'A' in vehicles and vehicles['A'][-1][1] == len(current_board[0]) - 1:
            return moves
        
        # Get blocking vehicles
        blocking = get_blocking_vehicles(current_board, vehicles)
        
        # Prioritize moving blocking vehicles
        for vehicle, positions in vehicles.items():
            # Skip non-blocking vehicles unless no blocking vehicles can move
            if vehicle != 'A' and vehicle not in blocking and blocking:
                continue
                
            for direction in [-1, 1]:  # Try moving up/left first
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}1"]
                        queue.append((new_board, new_moves))
                        
                        # If this move clears a direct path for AA, prioritize moving AA
                        if vehicle in blocking:
                            new_vehicles = get_vehicle_positions(new_board)
                            if not get_blocking_vehicles(new_board, new_vehicles):
                                if can_move(new_board, new_vehicles['A'], 1):
                                    final_board = move_vehicle(new_board, 'A', new_vehicles['A'], 1)
                                    return new_moves + ['A+1']
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")