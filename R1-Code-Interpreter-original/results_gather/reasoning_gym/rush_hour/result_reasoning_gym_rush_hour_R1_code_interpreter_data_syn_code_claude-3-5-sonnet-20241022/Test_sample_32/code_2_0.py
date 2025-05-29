from heapq import heappush, heappop
from collections import defaultdict

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board, vehicles):
    if 'A' not in vehicles:
        return float('inf')
    
    # Position of red car
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(x[1] for x in red_car)
    
    # Count blocking vehicles
    blocking = 0
    for vehicle, positions in vehicles.items():
        if vehicle != 'A':
            for pos in positions:
                if pos[0] == red_row and pos[1] > red_right:
                    blocking += 1
    
    return blocking

def get_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        is_horizontal = positions[0][0] == positions[1][0]
        if is_horizontal:
            row = positions[0][0]
            min_col = min(p[1] for p in positions)
            max_col = max(p[1] for p in positions)
            
            # Try left
            if min_col > 0 and board[row][min_col-1] == '.':
                moves.append((vehicle, -1))
            # Try right
            if max_col < len(board[0])-1 and board[row][max_col+1] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            min_row = min(p[0] for p in positions)
            max_row = max(p[0] for p in positions)
            
            # Try up
            if min_row > 0 and board[min_row-1][col] == '.':
                moves.append((vehicle, -1))
            # Try down
            if max_row < len(board)-1 and board[max_row+1][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = vehicles[vehicle]
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
        
    # Set new positions
    is_horizontal = positions[0][0] == positions[1][0]
    if is_horizontal:
        row = positions[0][0]
        cols = sorted(pos[1] for pos in positions)
        for col in cols:
            new_board[row][col + direction] = vehicle
    else:
        col = positions[0][1]
        rows = sorted(pos[0] for pos in positions)
        for row in rows:
            new_board[row + direction][col] = vehicle
            
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board):
    board = parse_board(initial_board)
    initial_vehicles = get_vehicle_positions(board)
    
    # Priority queue for A* search
    queue = [(0, 0, board, [])]  # (priority, moves_count, board, moves)
    seen = {board_to_string(board)}
    
    while queue:
        _, moves_count, current_board, moves = heappop(queue)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if 'A' in vehicles:
            red_car = vehicles['A']
            if max(pos[1] for pos in red_car) == len(current_board[0])-1:
                return moves
        
        # Get possible moves
        possible_moves = get_moves(current_board, vehicles)
        for vehicle, direction in possible_moves:
            new_board = apply_move(current_board, vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_vehicles = get_vehicle_positions(new_board)
                h_score = heuristic(new_board, new_vehicles)
                priority = moves_count + 1 + h_score
                
                move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                heappush(queue, (priority, moves_count + 1, new_board, moves + [move_str]))
    
    return None

board = """...HI.
BBBHI.
AA.HJ.
F.CCJ.
F.GDD.
EEG..."""

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")