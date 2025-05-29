from heapq import heappush, heappop

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

def blocking_vehicles(board, vehicles):
    if 'A' not in vehicles:
        return float('inf')
    
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(x[1] for x in red_car)
    
    blockers = set()
    for col in range(red_right + 1, len(board[0])):
        if board[red_row][col] != '.' and board[red_row][col] != 'x':
            blockers.add(board[red_row][col])
            
    # Add vehicles blocking the blockers
    secondary_blockers = set()
    for blocker in blockers:
        pos = vehicles[blocker]
        if pos[0][0] == pos[1][0]:  # horizontal
            continue
        col = pos[0][1]
        min_row = min(p[0] for p in pos)
        max_row = max(p[0] for p in pos)
        for row in range(min_row, max_row + 1):
            if board[row][col] != blocker and board[row][col] != '.':
                secondary_blockers.add(board[row][col])
    
    return len(blockers) + len(secondary_blockers) * 0.5

def get_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        is_horizontal = positions[0][0] == positions[1][0]
        if is_horizontal:
            row = positions[0][0]
            min_col = min(p[1] for p in positions)
            max_col = max(p[1] for p in positions)
            
            if min_col > 0 and board[row][min_col-1] == '.':
                moves.append((vehicle, -1))
            if max_col < len(board[0])-1 and board[row][max_col+1] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            min_row = min(p[0] for p in positions)
            max_row = max(p[0] for p in positions)
            
            if min_row > 0 and board[min_row-1][col] == '.':
                moves.append((vehicle, -1))
            if max_row < len(board)-1 and board[max_row+1][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
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
        
        possible_moves = get_moves(current_board, vehicles)
        for vehicle, direction in possible_moves:
            new_board = apply_move(current_board, vehicle, vehicles[vehicle], direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_vehicles = get_vehicle_positions(new_board)
                h_score = blocking_vehicles(new_board, new_vehicles)
                
                # Prioritize moves of blocking vehicles and the red car
                priority = moves_count + 1 + h_score
                if vehicle == 'A' or vehicle in [board[vehicles['A'][0][0]][col] 
                    for col in range(max(p[1] for p in vehicles['A'])+1, len(board[0]))
                    if board[vehicles['A'][0][0]][col] != '.']:
                    priority -= 0.5
                
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