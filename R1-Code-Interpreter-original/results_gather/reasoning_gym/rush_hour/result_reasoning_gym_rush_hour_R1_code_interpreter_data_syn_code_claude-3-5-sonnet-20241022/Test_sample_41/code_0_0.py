from collections import defaultdict, deque
import copy

def parse_board(board_str):
    # Convert board string to 2D list
    board = [list(row) for row in board_str.strip().split('\n')]
    
    # Find all vehicles and their positions
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                vehicles[board[i][j]].append((i, j))
    
    # Determine orientation for each vehicle
    orientations = {}
    for vehicle, positions in vehicles.items():
        if positions[0][0] == positions[-1][0]:  # same row
            orientations[vehicle] = 'H'  # horizontal
        else:
            orientations[vehicle] = 'V'  # vertical
    
    return board, vehicles, orientations

def get_valid_moves(board, vehicles, orientations):
    moves = []
    for vehicle in vehicles:
        positions = vehicles[vehicle]
        orientation = orientations[vehicle]
        
        if orientation == 'H':
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((vehicle, -1))
            
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((vehicle, -1))
            
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    new_vehicles = copy.deepcopy(vehicles)
    
    # Clear old positions
    for pos in vehicles[vehicle]:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    if vehicle in orientations and orientations[vehicle] == 'H':
        new_positions = [(p[0], p[1] + direction) for p in vehicles[vehicle]]
    else:
        new_positions = [(p[0] + direction, p[1]) for p in vehicles[vehicle]]
    
    # Update new positions
    new_vehicles[vehicle] = new_positions
    for pos in new_positions:
        new_board[pos[0]][pos[1]] = vehicle
    
    return new_board, new_vehicles

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board, vehicles, orientations = parse_board(initial_board_str)
    start_state = (board_to_string(board), vehicles)
    visited = {start_state[0]}
    queue = deque([(start_state, [])])
    
    while queue:
        (board_str, current_vehicles), moves = queue.popleft()
        board = [list(row) for row in board_str.split('\n')]
        
        # Check if red car (AA) reached the exit
        aa_positions = current_vehicles['A']
        if max(p[1] for p in aa_positions) == len(board[0])-1:
            return moves
        
        # Try all possible moves
        valid_moves = get_valid_moves(board, current_vehicles, orientations)
        for vehicle, direction in valid_moves:
            new_board, new_vehicles = apply_move(board, current_vehicles, vehicle, direction)
            new_board_str = board_to_string(new_board)
            
            if new_board_str not in visited:
                visited.add(new_board_str)
                new_moves = moves + [(vehicle, direction)]
                queue.append(((new_board_str, new_vehicles), new_moves))
    
    return None

# Initial board
initial_board = """x.IBB.
..IJCC
AAIJK.
.HDDKL
.HEEEL
FFGG.."""

# Solve the puzzle
solution = solve_puzzle(initial_board)

# Format the solution
if solution:
    formatted_solution = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution.append(f"{vehicle}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_solution) + ">>>")
else:
    print("<<<No solution found>>>")