from collections import defaultdict, deque
import copy

def parse_board(board_str):
    board = [list(row) for row in board_str.split('\n') if row]
    return board

def find_vehicles(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                vehicles[board[i][j]].append((i, j))
    return vehicles

def get_orientation(coords):
    if coords[0][0] == coords[1][0]:  # same row
        return 'horizontal'
    return 'vertical'

def is_valid_move(board, vehicle_coords, direction, steps):
    orientation = get_orientation(vehicle_coords)
    new_coords = []
    
    for x, y in vehicle_coords:
        if orientation == 'horizontal':
            new_y = y + steps
            if new_y < 0 or new_y >= len(board[0]) or board[x][new_y] not in ['.', vehicle_coords[0][1]]:
                return False
            new_coords.append((x, new_y))
        else:  # vertical
            new_x = x + steps
            if new_x < 0 or new_x >= len(board) or board[new_x][y] not in ['.', vehicle_coords[0][0]]:
                return False
            new_coords.append((new_x, y))
    return True

def make_move(board, vehicle, coords, steps):
    new_board = [row[:] for row in board]
    orientation = get_orientation(coords)
    
    # Clear old position
    for x, y in coords:
        new_board[x][y] = '.'
    
    # Set new position
    for x, y in coords:
        if orientation == 'horizontal':
            new_board[x][y + steps] = vehicle
        else:
            new_board[x + steps][y] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = parse_board("BBH...\nx.HCCC\n.AAI.K\n.G.IJK\n.GDDJK\n.EEEFF")
    vehicles = find_vehicles(initial_board)
    
    # BFS
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) reached the exit
        vehicles = find_vehicles(current_board)
        red_car = vehicles['A']
        if red_car[1][1] == len(current_board[0]) - 1:
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, coords in vehicles.items():
            orientation = get_orientation(coords)
            
            # Try moving in both directions
            for steps in [-1, 1]:
                if is_valid_move(current_board, coords, orientation, steps):
                    new_board = make_move(current_board, vehicle, coords, steps)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                        queue.append((new_board, new_moves))
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")