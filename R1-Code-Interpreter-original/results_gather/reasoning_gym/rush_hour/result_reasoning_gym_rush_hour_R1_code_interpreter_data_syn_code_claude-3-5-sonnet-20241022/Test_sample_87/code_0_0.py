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
        for pos in positions:
            new_board[pos[0]][pos[1] + direction] = vehicle
    else:
        for pos in positions:
            new_board[pos[0] + direction][pos[1]] = vehicle
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = parse_board("""G..IBB
G.HI.J
G.HAAJ
CCDD.K
EEEFFK
......""")
    
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car 'AA' reaches right edge)
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles['A']):
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:  # Try moving both directions
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")