from collections import deque, defaultdict

def parse_board(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_blocking_vehicles(board, vehicles):
    # Find vehicles blocking the red car (A)
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(x[1] for x in red_car)
    blocking = set()
    
    # Check all positions to the right of the red car
    for j in range(red_right + 1, len(board[0])):
        if board[red_row][j] not in '.x':
            blocking.add(board[red_row][j])
    return blocking

def get_moves(board, vehicles):
    moves = []
    height, width = len(board), len(board[0])
    
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try left
            if left > 0 and board[row][left-1] == '.':
                moves.append((vehicle, -1))
            # Try right
            if right < width-1 and board[row][right+1] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Try up
            if top > 0 and board[top-1][col] == '.':
                moves.append((vehicle, -1))
            # Try down
            if bottom < height-1 and board[bottom+1][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, move):
    vehicle, direction = move
    new_board = [list(row) for row in board]
    positions = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == vehicle]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    for i, j in positions:
        new_i = i + (0 if horizontal else direction)
        new_j = j + (direction if horizontal else 0)
        new_board[new_i][new_j] = vehicle
    
    return [''.join(row) for row in new_board]

def solve_rush_hour(board):
    initial_state = tuple(board)
    queue = deque([(initial_state, [])])
    seen = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()
        current_board = list(current_state)
        vehicles = parse_board(current_board)
        
        # Check if solved
        red_car = vehicles['A']
        if max(x[1] for x in red_car) == len(board[0])-1:
            return path
        
        # Get blocking vehicles and prioritize their moves
        blocking = get_blocking_vehicles(current_board, vehicles)
        moves = get_moves(current_board, vehicles)
        
        # Prioritize moves of blocking vehicles
        moves.sort(key=lambda m: m[0] in blocking, reverse=True)
        
        for move in moves:
            new_board = apply_move(current_board, move)
            new_state = tuple(new_board)
            
            if new_state not in seen:
                seen.add(new_state)
                move_str = f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"
                queue.append((new_state, path + [move_str]))
    
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