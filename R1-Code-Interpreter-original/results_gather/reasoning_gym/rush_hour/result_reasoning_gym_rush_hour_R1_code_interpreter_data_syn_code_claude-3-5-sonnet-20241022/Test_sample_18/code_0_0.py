from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, vehicle_pos, direction):
    if is_horizontal(vehicle_pos):
        if direction > 0:  # Move right
            new_pos = vehicle_pos[-1][1] + 1
            return new_pos < len(board[0]) and board[vehicle_pos[0][0]][new_pos] == '.'
        else:  # Move left
            new_pos = vehicle_pos[0][1] - 1
            return new_pos >= 0 and board[vehicle_pos[0][0]][new_pos] == '.'
    else:
        if direction > 0:  # Move down
            new_pos = vehicle_pos[-1][0] + 1
            return new_pos < len(board) and board[new_pos][vehicle_pos[0][1]] == '.'
        else:  # Move up
            new_pos = vehicle_pos[0][0] - 1
            return new_pos >= 0 and board[new_pos][vehicle_pos[0][1]] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        start_col = positions[0][1] + direction
        for i in range(len(positions)):
            new_board[row][start_col + i] = vehicle
    else:
        col = positions[0][1]
        start_row = positions[0][0] + direction
        for i in range(len(positions)):
            new_board[start_row + i][col] = vehicle
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    queue = deque([(board, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (AA reaches right edge)
        if 'A' in vehicles and vehicles['A'][-1][1] == len(current_board[0]) - 2:
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    return None

initial_board = """BBHICC
.GHIJ.
.GAAJK
DD..xK
EE...K
..xFF."""

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")