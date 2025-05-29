from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def find_vehicle_coords(board, vehicle):
    coords = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == vehicle:
                coords.append((i, j))
    return sorted(coords)

def is_horizontal(coords):
    return coords[0][0] == coords[1][0]

def can_move(board, coords, delta_row, delta_col):
    if not (0 <= coords[0][0] + delta_row < len(board) and 
            0 <= coords[-1][0] + delta_row < len(board)):
        return False
    
    if not (0 <= coords[0][1] + delta_col < len(board[0]) and 
            0 <= coords[-1][1] + delta_col < len(board[0])):
        return False
    
    vehicle = board[coords[0][0]][coords[0][1]]
    for coord in coords:
        new_row, new_col = coord[0] + delta_row, coord[1] + delta_col
        if board[new_row][new_col] not in ['.', vehicle]:
            return False
    return True

def make_move(board, coords, delta_row, delta_col):
    new_board = [row[:] for row in board]
    vehicle = board[coords[0][0]][coords[0][1]]
    
    # Clear old positions
    for row, col in coords:
        new_board[row][col] = '.'
    
    # Set new positions
    for row, col in coords:
        new_board[row + delta_row][col + delta_col] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def get_possible_moves(board, vehicle):
    coords = find_vehicle_coords(board, vehicle)
    moves = []
    
    if is_horizontal(coords):
        # Try left and right moves
        for steps in range(-3, 4):
            if steps == 0:
                continue
            if can_move(board, coords, 0, steps):
                moves.append((0, steps))
    else:
        # Try up and down moves
        for steps in range(-3, 4):
            if steps == 0:
                continue
            if can_move(board, coords, steps, 0):
                moves.append((steps, 0))
    
    return moves

def solve_puzzle():
    initial_board = parse_board("BBH...\nx.HCCC\n.AAI.K\n.G.IJK\n.GDDJK\n.EEEFF")
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    vehicles = set(c for row in initial_board for c in row if c not in ['.', 'x'])
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        red_car_coords = find_vehicle_coords(current_board, 'A')
        if max(col for _, col in red_car_coords) == len(current_board[0]) - 1:
            return moves
        
        # Try moves for each vehicle
        for vehicle in vehicles:
            coords = find_vehicle_coords(current_board, vehicle)
            possible_moves = get_possible_moves(current_board, vehicle)
            
            for delta_row, delta_col in possible_moves:
                new_board = make_move(current_board, coords, delta_row, delta_col)
                board_str = board_to_string(new_board)
                
                if board_str not in seen:
                    seen.add(board_str)
                    move_str = f"{vehicle}{'+' if delta_col > 0 or delta_row > 0 else ''}{delta_col if delta_col != 0 else delta_row}"
                    new_moves = moves + [move_str]
                    queue.append((new_board, new_moves))
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")