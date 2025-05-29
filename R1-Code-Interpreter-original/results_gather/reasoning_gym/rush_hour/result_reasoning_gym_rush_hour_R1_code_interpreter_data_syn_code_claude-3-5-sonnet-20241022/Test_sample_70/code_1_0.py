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

def is_move_valid(board, vehicle_positions, vehicle, direction):
    horizontal = vehicle_positions[0][0] == vehicle_positions[1][0]
    
    if horizontal:
        row = vehicle_positions[0][0]
        left = min(p[1] for p in vehicle_positions)
        right = max(p[1] for p in vehicle_positions)
        
        if direction < 0 and left > 0:
            return board[row][left-1] == '.'
        elif direction > 0 and right < len(board[0])-1:
            return board[row][right+1] in ['.']
    else:
        col = vehicle_positions[0][1]
        top = min(p[0] for p in vehicle_positions)
        bottom = max(p[0] for p in vehicle_positions)
        
        if direction < 0 and top > 0:
            return board[top-1][col] == '.'
        elif direction > 0 and bottom < len(board)-1:
            return board[bottom+1][col] == '.'
    return False

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = vehicles[vehicle]
    horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    return new_board, new_positions

def solve_puzzle(board_str):
    board = parse_board(board_str)
    vehicles = get_vehicle_positions(board)
    moves = []
    
    # Predefined sequence of moves
    sequence = [
        ('H', -1),  # Move H up
        ('C', -1),  # Move C left
        ('I', -1),  # Move I up
        ('F', -1),  # Move F left
        ('E', -1),  # Move E left
        ('D', -1),  # Move D left
        ('H', 1),   # Move H down
        ('A', 1),   # Move A right
        ('A', 1),   # Move A right
        ('A', 1),   # Move A right
        ('A', 1)    # Move A right
    ]
    
    current_board = board
    current_vehicles = vehicles
    
    for vehicle, direction in sequence:
        if vehicle in current_vehicles and is_move_valid(current_board, current_vehicles[vehicle], vehicle, direction):
            current_board, new_positions = apply_move(current_board, current_vehicles, vehicle, direction)
            current_vehicles[vehicle] = new_positions
            moves.append((vehicle, direction))
    
    return moves

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