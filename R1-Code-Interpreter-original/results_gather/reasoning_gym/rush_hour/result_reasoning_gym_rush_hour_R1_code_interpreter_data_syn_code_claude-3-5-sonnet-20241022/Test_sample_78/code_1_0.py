from collections import deque

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    # Sort positions for consistency
    for v in vehicles:
        vehicles[v].sort()
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def make_move(board, vehicle, positions, delta):
    new_board = [row[:] for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in positions:
        if horizontal:
            new_board[i][j + delta] = vehicle
        else:
            new_board[i + delta][j] = vehicle
            
    return new_board

def can_move(board, positions, delta):
    horizontal = is_horizontal(positions)
    if horizontal:
        row = positions[0][0]
        if delta > 0:  # Moving right
            rightmost = max(p[1] for p in positions)
            return (rightmost + delta < 6 and 
                   all(0 <= p[1] + delta < 6 and board[row][p[1] + delta] == '.' 
                       for p in positions))
        else:  # Moving left
            leftmost = min(p[1] for p in positions)
            return (leftmost + delta >= 0 and 
                   all(0 <= p[1] + delta < 6 and board[row][p[1] + delta] == '.' 
                       for p in positions))
    else:  # Vertical
        col = positions[0][1]
        if delta > 0:  # Moving down
            bottommost = max(p[0] for p in positions)
            return (bottommost + delta < 6 and 
                   all(0 <= p[0] + delta < 6 and board[p[0] + delta][col] == '.' 
                       for p in positions))
        else:  # Moving up
            topmost = min(p[0] for p in positions)
            return (topmost + delta >= 0 and 
                   all(0 <= p[0] + delta < 6 and board[p[0] + delta][col] == '.' 
                       for p in positions))

def is_solved(board):
    return board[2][4:6] == ['A', 'A']

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_blocking_vehicles(board, vehicles):
    # Find vehicles blocking the red car (AA)
    aa_pos = vehicles['A']
    aa_row = aa_pos[0][0]
    aa_rightmost = max(p[1] for p in aa_pos)
    
    blocking = set()
    for col in range(aa_rightmost + 1, 6):
        if board[aa_row][col] not in '.x':
            blocking.add(board[aa_row][col])
    return blocking

def solve_puzzle():
    initial_board = [
        list("BBBH.."),
        list("F..H.."),
        list("FAAI.J"),
        list("CCGI.J"),
        list("..GDDJ"),
        list("xxEE..")
    ]
    
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        blocking = get_blocking_vehicles(current_board, vehicles)
        
        # Prioritize moving blocking vehicles and the red car
        vehicle_order = list(blocking) + ['A'] + [v for v in vehicles if v not in blocking and v != 'A']
        
        for vehicle in vehicle_order:
            positions = vehicles[vehicle]
            is_horiz = is_horizontal(positions)
            
            # Try smaller moves first
            for delta in [1, -1, 2, -2, 3, -3]:
                if can_move(current_board, positions, delta):
                    new_board = make_move(current_board, vehicle, positions, delta)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if delta > 0 else '-'}{abs(delta)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")