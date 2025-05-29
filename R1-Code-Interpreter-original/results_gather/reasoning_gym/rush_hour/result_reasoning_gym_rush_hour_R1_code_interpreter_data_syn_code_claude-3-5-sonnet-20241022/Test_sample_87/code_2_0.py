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

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def calculate_heuristic(board, vehicles):
    # Calculate minimum moves needed to get red car to exit
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(pos[1] for pos in red_car)
    blocking_count = 0
    distance_to_exit = len(board[0]) - red_right - 1
    
    # Count blocking vehicles and their dependencies
    for col in range(red_right + 1, len(board[0])):
        if board[red_row][col] != '.':
            blocking_count += 2  # Each blocking vehicle needs at least 2 moves
    
    return blocking_count + distance_to_exit

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    is_horiz = is_horizontal(positions)
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    for row, col in positions:
        new_row = row + (0 if is_horiz else direction)
        new_col = col + (direction if is_horiz else 0)
        new_board[new_row][new_col] = vehicle
    
    return new_board

def can_move(board, positions, direction):
    is_horiz = is_horizontal(positions)
    if is_horiz:
        row = positions[0][0]
        if direction > 0:
            col = max(pos[1] for pos in positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:
            col = min(pos[1] for pos in positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = positions[0][1]
        if direction > 0:
            row = max(pos[0] for pos in positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:
            row = min(pos[0] for pos in positions) - 1
            return row >= 0 and board[row][col] == '.'

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = parse_board("""G..IBB
G.HI.J
G.HAAJ
CCDD.K
EEEFFK
......""")
    
    initial_vehicles = get_vehicle_positions(initial_board)
    initial_h = calculate_heuristic(initial_board, initial_vehicles)
    
    # Priority queue: (f_score, moves_count, board_string, moves)
    pq = [(initial_h, 0, board_to_string(initial_board), initial_board, [])]
    seen = {board_to_string(initial_board)}
    
    while pq:
        _, moves_count, _, current_board, moves = heappop(pq)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles['A']):
            return moves
        
        # Try all possible moves
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        new_vehicles = get_vehicle_positions(new_board)
                        h = calculate_heuristic(new_board, new_vehicles)
                        f_score = len(new_moves) + h
                        
                        heappush(pq, (f_score, len(new_moves), board_str, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")