from collections import deque

def get_board_state(board):
    # Convert board to a more manageable state
    state = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in state:
                    state[board[i][j]] = [(i, j)]
                else:
                    state[board[i][j]].append((i, j))
    return state

def is_horizontal(coords):
    return coords[0][0] == coords[1][0]

def can_move_right(vehicle_coords, board):
    row = vehicle_coords[0][0]
    rightmost = max(x[1] for x in vehicle_coords)
    return rightmost + 1 < len(board[0]) and board[row][rightmost + 1] == '.'

def can_move_left(vehicle_coords, board):
    row = vehicle_coords[0][0]
    leftmost = min(x[1] for x in vehicle_coords)
    return leftmost > 0 and board[row][leftmost - 1] == '.'

def can_move_up(vehicle_coords, board):
    col = vehicle_coords[0][1]
    topmost = min(x[0] for x in vehicle_coords)
    return topmost > 0 and board[topmost - 1][col] == '.'

def can_move_down(vehicle_coords, board):
    col = vehicle_coords[0][1]
    bottommost = max(x[0] for x in vehicle_coords)
    return bottommost + 1 < len(board) and board[bottommost + 1][col] == '.'

def move_vehicle(board, vehicle, direction):
    new_board = [list(row) for row in board]
    coords = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == vehicle]
    
    # Clear current position
    for i, j in coords:
        new_board[i][j] = '.'
    
    # Move to new position
    for i, j in coords:
        if direction == 'right':
            new_board[i][j + 1] = vehicle
        elif direction == 'left':
            new_board[i][j - 1] = vehicle
        elif direction == 'up':
            new_board[i - 1][j] = vehicle
        else:  # down
            new_board[i + 1][j] = vehicle
    
    return [''.join(row) for row in new_board]

def solve_rush_hour():
    initial_board = [
        "BBH.CC",
        "G.H.JK",
        "G.AAJK",
        "DD.IxL",
        "EE.I.L",
        "FFF..x"
    ]
    
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        state = get_board_state(current_board)
        
        # Check if solved (AA reaches exit)
        if 'A' in state and max(coord[1] for coord in state['A']) == 5:
            return moves
        
        # Focus on vehicles blocking AA's path first
        red_car_row = state['A'][0][0]
        red_car_right = max(coord[1] for coord in state['A'])
        blocking_vehicles = []
        
        for vehicle, coords in state.items():
            if vehicle != 'A':
                if is_horizontal(coords):
                    if coords[0][0] == red_car_row and min(coord[1] for coord in coords) > red_car_right:
                        blocking_vehicles.append(vehicle)
                else:  # vertical
                    if any(coord[0] == red_car_row and coord[1] > red_car_right for coord in coords):
                        blocking_vehicles.append(vehicle)
        
        # Prioritize moving blocking vehicles
        for vehicle in blocking_vehicles + list(state.keys()):
            coords = state[vehicle]
            if is_horizontal(coords):
                if can_move_right(coords, current_board):
                    new_board = move_vehicle(current_board, vehicle, 'right')
                    if tuple(new_board) not in seen:
                        seen.add(tuple(new_board))
                        queue.append((new_board, moves + [(vehicle, 1)]))
                if can_move_left(coords, current_board):
                    new_board = move_vehicle(current_board, vehicle, 'left')
                    if tuple(new_board) not in seen:
                        seen.add(tuple(new_board))
                        queue.append((new_board, moves + [(vehicle, -1)]))
            else:
                if can_move_up(coords, current_board):
                    new_board = move_vehicle(current_board, vehicle, 'up')
                    if tuple(new_board) not in seen:
                        seen.add(tuple(new_board))
                        queue.append((new_board, moves + [(vehicle, -1)]))
                if can_move_down(coords, current_board):
                    new_board = move_vehicle(current_board, vehicle, 'down')
                    if tuple(new_board) not in seen:
                        seen.add(tuple(new_board))
                        queue.append((new_board, moves + [(vehicle, 1)]))

solution = solve_rush_hour()
if solution:
    formatted_moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")