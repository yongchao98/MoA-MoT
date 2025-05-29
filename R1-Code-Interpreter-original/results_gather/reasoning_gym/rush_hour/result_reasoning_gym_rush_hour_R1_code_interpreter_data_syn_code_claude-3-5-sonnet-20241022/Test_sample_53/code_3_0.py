from collections import deque
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

def is_blocking_red_car(vehicle_pos, red_car_pos):
    if not vehicle_pos or not red_car_pos:
        return False
    return (vehicle_pos[0][0] == red_car_pos[0][0] and  # same row
            vehicle_pos[0][1] > max(pos[1] for pos in red_car_pos))  # to the right of red car

def get_blocking_vehicles(vehicles):
    red_car_pos = vehicles.get('A', [])
    blocking = []
    for vehicle, positions in vehicles.items():
        if vehicle != 'A' and is_blocking_red_car(positions, red_car_pos):
            blocking.append(vehicle)
    return blocking

def can_move(board, vehicle_positions, direction):
    if len(vehicle_positions) < 2:
        return False
    is_horiz = vehicle_positions[0][0] == vehicle_positions[1][0]
    if is_horiz:
        row = vehicle_positions[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in vehicle_positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in vehicle_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = vehicle_positions[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in vehicle_positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in vehicle_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [list(row) for row in board]
    is_horiz = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if is_horiz:
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = vehicle
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = vehicle
    
    return new_board

def solve_puzzle():
    initial_board = [
        list("..Ex.."),
        list("..EG.I"),
        list("AAFGHI"),
        list(".DF.H."),
        list(".DBBH."),
        list("..CCC.")
    ]
    
    queue = deque([(initial_board, [])])
    seen = {get_board_string(initial_board)}
    
    # Predefined sequence of priorities
    priority_sequence = ['F', 'H', 'B', 'C', 'A']
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles.get('A', [])):
            return moves
        
        # Try moves in priority order
        for priority_vehicle in priority_sequence:
            if priority_vehicle in vehicles:
                positions = vehicles[priority_vehicle]
                for direction in [1, -1]:
                    if can_move(current_board, positions, direction):
                        new_board = move_vehicle(current_board, priority_vehicle, positions, direction)
                        board_string = get_board_string(new_board)
                        
                        if board_string not in seen:
                            seen.add(board_string)
                            new_moves = moves + [f"{priority_vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                            
                            # Prioritize moves that follow the known solution pattern
                            if len(new_moves) <= 5:  # Only consider short solutions
                                queue.appendleft((new_board, new_moves))
                            else:
                                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")