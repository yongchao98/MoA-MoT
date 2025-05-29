from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def can_move(board, vehicle_positions, direction):
    if is_horizontal(vehicle_positions):
        row = vehicle_positions[0][0]
        if direction > 0:  # moving right
            col = max(pos[1] for pos in vehicle_positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # moving left
            col = min(pos[1] for pos in vehicle_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = vehicle_positions[0][1]
        if direction > 0:  # moving down
            row = max(pos[0] for pos in vehicle_positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # moving up
            row = min(pos[0] for pos in vehicle_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        cols = [pos[1] for pos in positions]
        new_cols = [col + direction for col in cols]
        for col in new_cols:
            new_board[row][col] = vehicle
    else:
        col = positions[0][1]
        rows = [pos[0] for pos in positions]
        new_rows = [row + direction for row in rows]
        for row in new_rows:
            new_board[row][col] = vehicle
    return new_board

def solve_puzzle(initial_board):
    initial_board = parse_board(initial_board)
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (red car AA reaches right edge)
        vehicles = get_vehicle_positions(current_board)
        if 'A' in vehicles:
            max_col = max(pos[1] for pos in vehicles['A'])
            if max_col == len(current_board[0]) - 1:
                return moves
        
        # Try all possible moves
        for vehicle, positions in vehicles.items():
            # Try moving in both directions
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = """...HI.
BBBHI.
AA.HJ.
F.CCJ.
F.GDD.
EEG..."""

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")