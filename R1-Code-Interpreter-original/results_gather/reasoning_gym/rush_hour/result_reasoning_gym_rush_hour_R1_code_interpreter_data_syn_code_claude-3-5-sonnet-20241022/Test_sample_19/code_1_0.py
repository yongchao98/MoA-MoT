from heapq import heappush, heappop

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_blocking_vehicles(board, vehicles, aa_row):
    blocking = set()
    if 'A' in vehicles:
        aa_right = max(p[1] for p in vehicles['A'])
        # Check direct horizontal blockers
        for col in range(aa_right + 1, len(board[0])):
            if board[aa_row][col] != '.' and board[aa_row][col] != 'x':
                blocking.add(board[aa_row][col])
    return blocking

def heuristic(board, vehicles, blocking):
    if 'A' not in vehicles:
        return float('inf')
    return len(blocking) * 2

def can_move(board, positions, direction):
    if is_horizontal(positions):
        if direction > 0:
            new_col = max(p[1] for p in positions) + 1
            return new_col < len(board[0]) and board[positions[0][0]][new_col] == '.'
        else:
            new_col = min(p[1] for p in positions) - 1
            return new_col >= 0 and board[positions[0][0]][new_col] == '.'
    else:
        if direction > 0:
            new_row = max(p[0] for p in positions) + 1
            return new_row < len(board) and board[new_row][positions[0][1]] == '.'
        else:
            new_row = min(p[0] for p in positions) - 1
            return new_row >= 0 and board[new_row][positions[0][1]] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    if is_horizontal(positions):
        for pos in positions:
            new_board[pos[0]][pos[1] + direction] = vehicle
    else:
        for pos in positions:
            new_board[pos[0] + direction][pos[1]] = vehicle
    return new_board

def solve_puzzle():
    initial_board = [
        ['.', 'F', 'B', 'B', 'J', 'K'],
        ['.', 'F', '.', 'I', 'J', 'K'],
        ['x', 'A', 'A', 'I', '.', 'L'],
        ['C', 'C', 'H', 'I', '.', 'L'],
        ['E', 'G', 'H', '.', '.', '.'],
        ['E', 'G', 'D', 'D', 'D', '.']
    ]
    
    aa_row = 2  # The row where AA is located
    visited = {board_to_string(initial_board)}
    pq = []  # Priority queue: (priority, moves_count, board, moves)
    
    vehicles = get_vehicle_positions(initial_board)
    blocking = get_blocking_vehicles(initial_board, vehicles, aa_row)
    heappush(pq, (heuristic(initial_board, vehicles, blocking), 0, initial_board, []))
    
    while pq:
        _, moves_count, current_board, moves = heappop(pq)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if 'A' in vehicles and max(p[1] for p in vehicles['A']) == len(current_board[0]) - 1:
            return moves
        
        blocking = get_blocking_vehicles(current_board, vehicles, aa_row)
        
        # First priority: Move blocking vehicles out of the way
        for vehicle, positions in vehicles.items():
            # Prioritize moving blocking vehicles
            is_priority = vehicle in blocking
            
            for direction in [-1, 1]:  # Try moving up/left first
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}1"]
                        new_vehicles = get_vehicle_positions(new_board)
                        new_blocking = get_blocking_vehicles(new_board, new_vehicles, aa_row)
                        
                        # Lower priority for non-blocking moves
                        priority = heuristic(new_board, new_vehicles, new_blocking)
                        if not is_priority:
                            priority += 10
                            
                        heappush(pq, (priority, moves_count + 1, new_board, new_moves))
                        
                        # If path is clear for AA, prioritize moving it
                        if vehicle in blocking and not new_blocking:
                            if can_move(new_board, new_vehicles['A'], 1):
                                return new_moves + ['A+1']
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")