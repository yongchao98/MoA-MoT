from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n')]

def get_vehicle_info(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = {'positions': [], 'orientation': None}
                vehicles[board[i][j]]['positions'].append((i, j))
    
    for v in vehicles:
        if vehicles[v]['positions'][0][0] == vehicles[v]['positions'][1][0]:
            vehicles[v]['orientation'] = 'horizontal'
        else:
            vehicles[v]['orientation'] = 'vertical'
    return vehicles

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def make_move(board, vehicle, positions, direction, steps):
    new_board = [row[:] for row in board]
    for i, j in positions:
        new_board[i][j] = '.'
    
    if direction == 'horizontal':
        new_positions = [(i, j + steps) for i, j in positions]
    else:
        new_positions = [(i + steps, j) for i, j in positions]
        
    for i, j in new_positions:
        new_board[i][j] = vehicle
    return new_board

def is_valid_move(board, positions, direction, steps):
    if direction == 'horizontal':
        for i, j in positions:
            new_j = j + steps
            if new_j < 0 or new_j >= len(board[0]) or board[i][new_j] not in ['.', positions[0][1]]:
                return False
    else:
        for i, j in positions:
            new_i = i + steps
            if new_i < 0 or new_i >= len(board) or board[new_i][j] not in ['.', positions[0][0]]:
                return False
    return True

def solve_puzzle(initial_board):
    board = parse_board(initial_board)
    vehicles = get_vehicle_info(board)
    start_state = board_to_string(board)
    
    queue = deque([(start_state, [])])
    visited = {start_state}
    
    while queue:
        current_state, moves = queue.popleft()
        current_board = parse_board(current_state)
        
        # Check if red car (AA) can reach exit
        red_car = vehicles['A']['positions']
        if any(j == len(current_board[0])-1 for i, j in red_car):
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, info in vehicles.items():
            positions = info['positions']
            direction = info['orientation']
            
            # Try moving in both directions
            for steps in range(-5, 6):
                if steps == 0:
                    continue
                    
                if is_valid_move(current_board, positions, direction, steps):
                    new_board = make_move(current_board, vehicle, positions, direction, steps)
                    new_state = board_to_string(new_board)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                        queue.append((new_state, new_moves))
    
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