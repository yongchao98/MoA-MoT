from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.strip().split('\n')]

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_possible_moves(board, vehicles):
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
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            # Try moving up
            if top >= 0 and board[top][col] == '.':
                moves.append((vehicle, -1))
            
            # Try moving down
            if bottom < len(board) and board[bottom][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    return new_board

def solve_puzzle(initial_board):
    initial_board = parse_board(initial_board)
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car 'AA' reaches right edge)
        red_car = vehicles['A']
        if max(pos[1] for pos in red_car) == len(current_board[0]) - 2:
            return moves
        
        # Try all possible moves
        for vehicle, direction in get_possible_moves(current_board, vehicles):
            new_board = apply_move(current_board, vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(vehicle, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
board = """HBBCCx
H.IDDD
AAIJ..
.EEJFF
.GGG.K
.....K"""

# Solve the puzzle
solution = solve_puzzle(board)

# Format the solution
if solution:
    formatted_moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}1")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")