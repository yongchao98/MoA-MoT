from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.strip().split('\n')]

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

def can_move(board, vehicle_pos, direction):
    if is_horizontal(vehicle_pos):
        if direction > 0:  # Move right
            new_pos = vehicle_pos[-1][1] + 1
            return (new_pos < len(board[0]) and 
                   board[vehicle_pos[0][0]][new_pos] == '.')
        else:  # Move left
            new_pos = vehicle_pos[0][1] - 1
            return (new_pos >= 0 and 
                   board[vehicle_pos[0][0]][new_pos] == '.')
    else:  # Vertical
        if direction > 0:  # Move down
            new_pos = vehicle_pos[-1][0] + 1
            return (new_pos < len(board) and 
                   board[new_pos][vehicle_pos[0][1]] == '.')
        else:  # Move up
            new_pos = vehicle_pos[0][0] - 1
            return (new_pos >= 0 and 
                   board[new_pos][vehicle_pos[0][1]] == '.')

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        start_col = positions[0][1] + direction
        for col in range(start_col, start_col + len(positions)):
            new_board[row][col] = vehicle
    else:
        col = positions[0][1]
        start_row = positions[0][0] + direction
        for row in range(start_row, start_row + len(positions)):
            new_board[row][col] = vehicle
    return new_board

def solve_puzzle(initial_board):
    board = parse_board(initial_board)
    start_state = board_to_string(board)
    queue = deque([(start_state, [])])
    visited = {start_state}
    
    while queue:
        current_state, moves = queue.popleft()
        current_board = parse_board(current_state)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        if any(j == len(current_board[0])-1 for i, j in vehicles['A']):
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:  # Try moving both directions
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    new_state = board_to_string(new_board)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_state, new_moves))
    return None

initial_board = """x.BBCC
DDx.KL
I.AAKL
IEEFFL
GGGJ..
HHHJ.."""

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")