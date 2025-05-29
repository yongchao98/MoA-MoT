from collections import defaultdict, deque
import copy

def parse_board(board_str):
    # Convert string board to 2D list
    return [list(row) for row in board_str.split('\n') if row]

def find_vehicles(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                vehicles[board[i][j]].append((i, j))
    return vehicles

def get_orientation(coords):
    # Return 'h' for horizontal, 'v' for vertical
    return 'h' if coords[0][0] == coords[1][0] else 'v'

def get_valid_moves(board, vehicles):
    moves = []
    for vehicle, coords in vehicles.items():
        orientation = get_orientation(coords)
        if orientation == 'h':
            # Try moving left
            if coords[0][1] > 0 and board[coords[0][0]][coords[0][1]-1] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            if coords[-1][1] < len(board[0])-1 and board[coords[0][0]][coords[-1][1]+1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            if coords[0][0] > 0 and board[coords[0][0]-1][coords[0][1]] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            if coords[-1][0] < len(board)-1 and board[coords[-1][0]+1][coords[-1][1]] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    new_board = [row[:] for row in board]
    new_vehicles = copy.deepcopy(vehicles)
    coords = new_vehicles[vehicle]
    orientation = get_orientation(coords)
    
    # Clear current position
    for x, y in coords:
        new_board[x][y] = '.'
    
    # Update coordinates
    if orientation == 'h':
        new_vehicles[vehicle] = [(x, y + direction) for x, y in coords]
    else:
        new_vehicles[vehicle] = [(x + direction, y) for x, y in coords]
    
    # Update board
    for x, y in new_vehicles[vehicle]:
        new_board[x][y] = vehicle
    
    return new_board, new_vehicles

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = find_vehicles(board)
    
    # BFS
    queue = deque([(board, vehicles, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if red car (AA) is at exit
        red_car = current_vehicles['A']
        if red_car[-1][1] == len(current_board[0])-1:
            return moves
        
        # Try all possible moves
        for move in get_valid_moves(current_board, current_vehicles):
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                vehicle, spaces = move
                move_str = f"{vehicle}{'+' if spaces > 0 else ''}{spaces}"
                queue.append((new_board, new_vehicles, moves + [move_str]))
    
    return None

# Initial board
board_str = """.BBCCJ
x....J
GAAI.K
GDDI.K
..HIEE
FFH..."""

solution = solve_puzzle(board_str)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")