from collections import deque
import copy

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            # Try moving left
            if left >= 0 and board[row][left] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            if right < 6 and board[row][right] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            # Try moving up
            if top >= 0 and board[top][col] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            if bottom < 6 and board[bottom][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [list(row) for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board]

def is_solved(board):
    # Check if AA is at the exit (right side)
    return 'A' in board[2][-2:]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        valid_moves = get_valid_moves(current_board, vehicles)
        
        for vehicle, direction in valid_moves:
            new_board = apply_move(current_board, vehicles, vehicle, direction)
            new_board_tuple = tuple(new_board)
            
            if new_board_tuple not in seen:
                seen.add(new_board_tuple)
                new_moves = moves + [(vehicle, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
board = [
    "..IBB.",
    ".GI.K.",
    ".GAAK.",
    "FGJCC.",
    "FHJ...",
    "FHDDEE"
]

# Solve the puzzle
solution = solve_puzzle(board)

# Format the solution
if solution:
    formatted_solution = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution.append(f"{vehicle}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_solution) + ">>>")
else:
    print("<<<No solution found>>>")