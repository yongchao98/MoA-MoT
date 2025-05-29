from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

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

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def make_move(board, vehicle, positions, delta):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    # Set new positions
    if is_horizontal(positions):
        for i, j in positions:
            new_board[i][j + delta] = vehicle
    else:
        for i, j in positions:
            new_board[i + delta][j] = vehicle
    return new_board

def get_possible_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost - 1] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0]) - 1 and board[positions[0][0]][rightmost + 1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost - 1][positions[0][1]] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board) - 1 and board[bottommost + 1][positions[0][1]] == '.':
                moves.append((vehicle, 1))
    return moves

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    start_state = board_to_string(board)
    queue = deque([(start_state, [])])
    visited = {start_state}
    
    while queue:
        current_state, moves = queue.popleft()
        current_board = parse_board(current_state)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car 'AA' reaches right edge)
        aa_positions = vehicles['A']
        if max(p[1] for p in aa_positions) == len(current_board[0]) - 2:
            return moves
        
        # Try all possible moves
        possible_moves = get_possible_moves(current_board, vehicles)
        for vehicle, delta in possible_moves:
            new_board = make_move(current_board, vehicle, vehicles[vehicle], delta)
            new_state = board_to_string(new_board)
            
            if new_state not in visited:
                visited.add(new_state)
                new_moves = moves + [(vehicle, delta)]
                queue.append((new_state, new_moves))
    
    return None

# Initial board
board = """G.JBB.
G.JK.L
HAAK.L
HICC.x
HIDD..
EEEFF."""

# Solve and format output
solution = solve_puzzle(board)
if solution:
    formatted_moves = []
    for vehicle, delta in solution:
        sign = '+' if delta > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}{abs(delta)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")