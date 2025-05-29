from heapq import heappush, heappop
import copy

def get_board_string(board):
    return '\n'.join(''.join(row) for row in board)

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

def can_move(board, vehicle_positions, direction):
    if is_horizontal(vehicle_positions):
        row = vehicle_positions[0][0]
        if direction > 0:  # moving right
            col = max(pos[1] for pos in vehicle_positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # moving left
            col = min(pos[1] for pos in vehicle_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = vehicle_positions[0][1]
        if direction > 0:  # moving down
            row = max(pos[0] for pos in vehicle_positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # moving up
            row = min(pos[0] for pos in vehicle_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [list(row) for row in board]
    for row, col in positions:
        new_board[row][col] = '.'
    if is_horizontal(positions):
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = vehicle
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = vehicle
    return new_board

def calculate_heuristic(board, vehicles):
    # Get red car (A) position
    red_car = vehicles['A']
    red_car_row = red_car[0][0]
    red_car_rightmost = max(pos[1] for pos in red_car)
    
    # Count blocking vehicles
    blocking_count = 0
    for vehicle, positions in vehicles.items():
        if vehicle != 'A':
            for pos in positions:
                if pos[0] == red_car_row and pos[1] > red_car_rightmost:
                    blocking_count += 1
                    break
    
    return blocking_count

def solve_puzzle():
    initial_board = [
        list("..Ex.."),
        list("..EG.I"),
        list("AAFGHI"),
        list(".DF.H."),
        list(".DBBH."),
        list("..CCC.")
    ]
    
    # Priority queue for A* search
    pq = [(0, 0, initial_board, [])]  # (f_score, move_count, board, moves)
    seen = {get_board_string(initial_board)}
    move_count = 0
    
    while pq:
        _, current_moves_count, current_board, moves = heappop(pq)
        
        vehicles = get_vehicle_positions(current_board)
        # Check if solved
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles['A']):
            return moves
            
        # Try moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_string = get_board_string(new_board)
                    
                    if board_string not in seen:
                        seen.add(board_string)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        new_vehicles = get_vehicle_positions(new_board)
                        h_score = calculate_heuristic(new_board, new_vehicles)
                        g_score = len(new_moves)
                        f_score = g_score + h_score * 2  # Weighting heuristic more heavily
                        
                        heappush(pq, (f_score, g_score, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")