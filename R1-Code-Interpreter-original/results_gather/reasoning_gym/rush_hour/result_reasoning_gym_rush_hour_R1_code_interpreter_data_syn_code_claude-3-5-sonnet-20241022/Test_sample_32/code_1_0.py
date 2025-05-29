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

def evaluate_state(board, vehicles):
    if 'A' not in vehicles:
        return float('inf')
    
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(x[1] for x in red_car)
    
    # Direct blockers
    direct_blockers = set()
    for col in range(red_right + 1, len(board[0])):
        if board[red_row][col] != '.':
            direct_blockers.add(board[red_row][col])
    
    # Calculate space available for each blocking vehicle to move
    blocking_score = 0
    for blocker in direct_blockers:
        pos = vehicles[blocker]
        if pos[0][0] == pos[1][0]:  # horizontal blocker
            blocking_score += 3
        else:  # vertical blocker
            blocking_score += 2
            # Check if the blocker is itself blocked
            col = pos[0][1]
            min_row = min(p[0] for p in pos)
            max_row = max(p[0] for p in pos)
            blocked_spaces = 0
            for row in range(min_row - 1, max_row + 2):
                if row < 0 or row >= len(board) or board[row][col] != '.':
                    blocked_spaces += 1
            blocking_score += blocked_spaces * 0.5
    
    return blocking_score

def get_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        is_horizontal = positions[0][0] == positions[1][0]
        if is_horizontal:
            row = positions[0][0]
            min_col = min(p[1] for p in positions)
            max_col = max(p[1] for p in positions)
            
            if min_col > 0 and board[row][min_col-1] == '.':
                moves.append((vehicle, -1, 'h'))
            if max_col < len(board[0])-1 and board[row][max_col+1] == '.':
                moves.append((vehicle, 1, 'h'))
        else:
            col = positions[0][1]
            min_row = min(p[0] for p in positions)
            max_row = max(p[0] for p in positions)
            
            if min_row > 0 and board[min_row-1][col] == '.':
                moves.append((vehicle, -1, 'v'))
            if max_row < len(board)-1 and board[max_row+1][col] == '.':
                moves.append((vehicle, 1, 'v'))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction, orientation = move
    positions = vehicles[vehicle]
    new_board = [row[:] for row in board]
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if orientation == 'h':
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

def solve_puzzle(initial_board):
    board = parse_board(initial_board)
    initial_vehicles = get_vehicle_positions(board)
    
    queue = [(0, board, [])]  # (priority, board, moves)
    seen = {str(board)}
    
    while queue:
        _, current_board, moves = heappop(queue)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if 'A' in vehicles and max(pos[1] for pos in vehicles['A']) == len(current_board[0])-1:
            return moves
        
        possible_moves = get_moves(current_board, vehicles)
        for move in possible_moves:
            vehicle, direction, _ = move
            new_board = apply_move(current_board, vehicles, move)
            board_str = str(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_vehicles = get_vehicle_positions(new_board)
                score = evaluate_state(new_board, new_vehicles)
                
                # Prioritize moves that directly help clear the path
                if vehicle == 'A':
                    score -= 1
                elif 'A' in vehicles:
                    red_row = vehicles['A'][0][0]
                    red_right = max(pos[1] for pos in vehicles['A'])
                    if any(pos[0] == red_row and pos[1] > red_right for pos in vehicles[vehicle]):
                        score -= 0.5
                
                move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                heappush(queue, (score + len(moves), new_board, moves + [move_str]))
    
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