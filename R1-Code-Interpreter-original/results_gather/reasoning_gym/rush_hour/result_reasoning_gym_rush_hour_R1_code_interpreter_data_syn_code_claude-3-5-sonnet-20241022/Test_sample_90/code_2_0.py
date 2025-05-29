from collections import deque

def create_initial_state():
    # Define the board as a list of strings for easier manipulation
    board = [
        "BBBKCC",
        "DDJK.L",
        "IJAAL",
        "IEE.xM",
        "FF...M",
        "GGHH.x"
    ]
    
    # Find all vehicles and their positions
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            car = board[i][j]
            if car not in '.x':
                if car not in vehicles:
                    vehicles[car] = [(i, j)]
                else:
                    vehicles[car].append((i, j))
    
    # Determine orientation and length for each vehicle
    vehicle_info = {}
    for car, positions in vehicles.items():
        positions.sort()  # Sort positions to ensure consistent order
        is_horizontal = positions[0][0] == positions[-1][0]
        vehicle_info[car] = {
            'positions': positions,
            'horizontal': is_horizontal,
            'length': len(positions)
        }
    
    return board, vehicle_info

def can_move(board, car_info, direction):
    positions = car_info['positions']
    is_horizontal = car_info['horizontal']
    
    if is_horizontal:
        if direction > 0:  # Moving right
            row = positions[0][0]
            new_col = positions[-1][1] + 1
            if new_col >= len(board[0]) or board[row][new_col] not in '.':
                return False
        else:  # Moving left
            row = positions[0][0]
            new_col = positions[0][1] - 1
            if new_col < 0 or board[row][new_col] not in '.':
                return False
    else:
        if direction > 0:  # Moving down
            col = positions[0][1]
            new_row = positions[-1][0] + 1
            if new_row >= len(board) or board[new_row][col] not in '.':
                return False
        else:  # Moving up
            col = positions[0][1]
            new_row = positions[0][0] - 1
            if new_row < 0 or board[new_row][col] not in '.':
                return False
    return True

def move_vehicle(board, car, car_info, direction):
    new_board = [list(row) for row in board]
    new_positions = []
    
    # Clear current positions
    for pos in car_info['positions']:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    for pos in car_info['positions']:
        if car_info['horizontal']:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle():
    board, vehicle_info = create_initial_state()
    start_state = (board, vehicle_info)
    
    queue = deque([(start_state, [])])
    seen = {''.join(board)}
    
    while queue:
        (current_board, current_vehicles), moves = queue.popleft()
        
        # Check if red car (AA) is at the exit
        red_car = current_vehicles['A']
        if red_car['horizontal'] and red_car['positions'][-1][1] == 4:
            return moves
        
        # Try moving each vehicle
        for car, info in current_vehicles.items():
            directions = [-1, 1]
            for d in directions:
                if can_move(current_board, info, d):
                    new_board, new_positions = move_vehicle(current_board, car, info, d)
                    board_str = ''.join(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_vehicles = dict(current_vehicles)
                        new_vehicles[car] = dict(info)
                        new_vehicles[car]['positions'] = new_positions
                        
                        move = f"{car}{'+' if d > 0 else '-'}1"
                        queue.append(((new_board, new_vehicles), moves + [move]))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("<<<No solution found>>>")