def get_blocking_vehicles(board, row, start_col, end_col):
    blocking = set()
    for col in range(start_col + 1, end_col + 1):
        if board[row][col] not in ['.', 'x']:
            blocking.add(board[row][col])
    return blocking

def can_move(board, vehicle_coords, direction):
    if direction == 'up':
        return all(coord[0] > 0 and board[coord[0]-1][coord[1]] == '.' for coord in vehicle_coords)
    elif direction == 'down':
        return all(coord[0] < 5 and board[coord[0]+1][coord[1]] == '.' for coord in vehicle_coords)
    elif direction == 'left':
        return all(coord[1] > 0 and board[coord[0]][coord[1]-1] == '.' for coord in vehicle_coords)
    else:  # right
        return all(coord[1] < 5 and board[coord[0]][coord[1]+1] == '.' for coord in vehicle_coords)

def get_vehicle_coords(board, vehicle):
    coords = []
    for i in range(6):
        for j in range(6):
            if board[i][j] == vehicle:
                coords.append((i, j))
    return coords

def move_vehicle(board, vehicle, direction):
    board = [list(row) for row in board]
    coords = get_vehicle_coords(board, vehicle)
    
    # Clear current positions
    for i, j in coords:
        board[i][j] = '.'
    
    # Move to new positions
    for i, j in coords:
        new_i, new_j = i, j
        if direction == 'up': new_i -= 1
        elif direction == 'down': new_i += 1
        elif direction == 'left': new_j -= 1
        else: new_j += 1
        board[new_i][new_j] = vehicle
    
    return [''.join(row) for row in board]

def solve_puzzle():
    initial_board = [
        "BBH.CC",
        "G.H.JK",
        "G.AAJK",
        "DD.IxL",
        "EE.I.L",
        "FFF..x"
    ]
    
    # Find red car (AA) position
    aa_coords = get_vehicle_coords(initial_board, 'A')
    aa_row = aa_coords[0][0]
    aa_col = max(coord[1] for coord in aa_coords)
    
    # Get blocking vehicles
    blocking = get_blocking_vehicles(initial_board, aa_row, aa_col, 5)
    
    # We know J and K are blocking, and they can only move up
    # Try the direct solution
    current_board = initial_board
    moves = []
    
    # Move K up
    k_coords = get_vehicle_coords(current_board, 'K')
    if can_move(current_board, k_coords, 'up'):
        current_board = move_vehicle(current_board, 'K', 'up')
        moves.append(('K', -1))
    
    # Move J up
    j_coords = get_vehicle_coords(current_board, 'J')
    if can_move(current_board, j_coords, 'up'):
        current_board = move_vehicle(current_board, 'J', 'up')
        moves.append(('J', -1))
    
    # Move A right
    a_coords = get_vehicle_coords(current_board, 'A')
    if can_move(current_board, a_coords, 'right'):
        current_board = move_vehicle(current_board, 'A', 'right')
        moves.append(('A', 3))  # Move all the way to exit
    
    # Format the solution
    formatted_moves = []
    for vehicle, direction in moves:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}{abs(direction)}")
    
    return formatted_moves

solution = solve_puzzle()
print("<<<" + " ".join(solution) + ">>>")