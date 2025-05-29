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

def get_possible_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((vehicle, -1))
            
            # Try moving right
            if right < len(board[0])-1 and board[row][right+1] not in ['x', 'x']:
                if board[row][right+1] == '.':
                    moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((vehicle, -1))
            
            # Try moving down
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
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

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = get_vehicle_positions(board)
    
    queue = deque([(board, vehicles, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if solved (red car 'AA' reaches right edge)
        aa_positions = current_vehicles['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        possible_moves = get_possible_moves(current_board, current_vehicles)
        
        for vehicle, direction in possible_moves:
            new_board, new_positions = apply_move(current_board, current_vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_vehicles = dict(current_vehicles)
                new_vehicles[vehicle] = new_positions
                new_moves = moves + [(vehicle, direction)]
                queue.append((new_board, new_vehicles, new_moves))
    
    return None

# Initial board
board = """GBBB..
G.HCC.
AAH...
..HDDx
EEEI..
.FFI.."""

solution = solve_puzzle(board)
if solution:
    # Format the solution
    formatted_solution = ''
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{vehicle}{sign}{abs(direction)} "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")