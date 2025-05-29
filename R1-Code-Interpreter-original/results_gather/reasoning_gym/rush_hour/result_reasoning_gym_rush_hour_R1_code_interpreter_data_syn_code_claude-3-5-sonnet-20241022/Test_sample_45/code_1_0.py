from collections import deque

def parse_board(board):
    # Convert board to a more manageable format
    return [list(row) for row in board]

def get_vehicle_info(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_solved(board):
    # Check if red car (AA) reaches the right edge
    for i in range(len(board)):
        if board[i][-1] == 'A':
            return True
    return False

def get_moves(board, vehicles):
    moves = []
    height, width = len(board), len(board[0])
    
    for vehicle, coords in vehicles.items():
        # Determine if vehicle is horizontal or vertical
        is_horizontal = coords[0][0] == coords[1][0]
        
        if is_horizontal:
            row = coords[0][0]
            min_col = min(c[1] for c in coords)
            max_col = max(c[1] for c in coords)
            
            # Try move left
            if min_col > 0 and board[row][min_col - 1] == '.':
                moves.append((vehicle, -1))
            # Try move right
            if max_col < width - 1 and board[row][max_col + 1] == '.':
                moves.append((vehicle, 1))
        else:
            col = coords[0][1]
            min_row = min(c[0] for c in coords)
            max_row = max(c[0] for c in coords)
            
            # Try move up
            if min_row > 0 and board[min_row - 1][col] == '.':
                moves.append((vehicle, -1))
            # Try move down
            if max_row < height - 1 and board[max_row + 1][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    coords = vehicles[vehicle]
    is_horizontal = coords[0][0] == coords[1][0]
    
    # Clear current positions
    for i, j in coords:
        new_board[i][j] = '.'
    
    # Add new positions
    new_coords = []
    for i, j in coords:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = vehicle
        new_coords.append((new_i, new_j))
    
    return new_board, new_coords

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_rush_hour(board):
    initial_board = parse_board(board)
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        vehicles = get_vehicle_info(current_board)
        possible_moves = get_moves(current_board, vehicles)
        
        for vehicle, direction in possible_moves:
            new_board, new_coords = apply_move(current_board, vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                move_str = f"{vehicle}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = [
    "BBH.CC",
    ".GH..K",
    ".GAAJK",
    "DD.IJL",
    "...IxL",
    ".EEFFL"
]

solution = solve_rush_hour(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")