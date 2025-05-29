from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n')]

def get_vehicles(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            char = board[i][j]
            if char not in '.x':
                if char not in vehicles:
                    vehicles[char] = {'positions': [], 'orientation': None}
                vehicles[char]['positions'].append((i, j))
    
    # Determine orientation
    for v in vehicles:
        pos = vehicles[v]['positions']
        vehicles[v]['orientation'] = 'horizontal' if pos[0][0] == pos[1][0] else 'vertical'
    return vehicles

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def make_move(board, vehicle_char, positions, orientation, steps):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    if orientation == 'horizontal':
        for i, j in positions:
            new_board[i][j + steps] = vehicle_char
    else:
        for i, j in positions:
            new_board[i + steps][j] = vehicle_char
    return new_board

def is_valid_move(board, positions, orientation, steps):
    height, width = len(board), len(board[0])
    if orientation == 'horizontal':
        for i, j in positions:
            new_j = j + steps
            if new_j < 0 or new_j >= width:
                return False
            if board[i][new_j] not in ['.', board[i][j]] and (i, new_j) not in positions:
                return False
    else:
        for i, j in positions:
            new_i = i + steps
            if new_i < 0 or new_i >= height:
                return False
            if board[new_i][j] not in ['.', board[i][j]] and (new_i, j) not in positions:
                return False
    return True

def get_new_positions(positions, orientation, steps):
    if orientation == 'horizontal':
        return [(i, j + steps) for i, j in positions]
    return [(i + steps, j) for i, j in positions]

def solve_puzzle(board_str):
    board = parse_board(board_str)
    vehicles = get_vehicles(board)
    start_state = board_to_string(board)
    
    queue = deque([(board, [])])
    visited = {start_state}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) can reach exit
        red_car_pos = [(i, j) for i in range(len(current_board)) 
                      for j in range(len(current_board[0])) if current_board[i][j] == 'A']
        if max(j for _, j in red_car_pos) == len(current_board[0]) - 1:
            return moves
        
        # Try moves for each vehicle
        for vehicle, info in vehicles.items():
            orientation = info['orientation']
            current_pos = [(i, j) for i in range(len(current_board)) 
                          for j in range(len(current_board[0])) if current_board[i][j] == vehicle]
            
            # Try moving each direction
            max_steps = 5  # Maximum steps to try
            for steps in range(-max_steps, max_steps + 1):
                if steps == 0:
                    continue
                
                if is_valid_move(current_board, current_pos, orientation, steps):
                    new_board = make_move(current_board, vehicle, current_pos, orientation, steps)
                    new_state = board_to_string(new_board)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                        queue.append((new_board, new_moves))
                        
                        # Early exit if solution is found
                        if vehicle == 'A' and max(j for _, j in get_new_positions(current_pos, orientation, steps)) == len(current_board[0]) - 1:
                            return new_moves
    
    return None

# Initial board
board = """BBBCCM
..JDDM
AAJK..
IEEKFF
IGGKLx
HH..L."""

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")