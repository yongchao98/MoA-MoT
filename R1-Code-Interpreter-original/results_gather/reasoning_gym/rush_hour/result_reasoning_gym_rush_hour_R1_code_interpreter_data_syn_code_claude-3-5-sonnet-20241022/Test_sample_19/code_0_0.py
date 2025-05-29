from collections import deque
import copy

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, vehicle_pos, direction):
    if is_horizontal(vehicle_pos):
        if direction > 0:  # Move right
            new_col = vehicle_pos[-1][1] + 1
            return (new_col < len(board[0]) and 
                   board[vehicle_pos[0][0]][new_col] == '.')
        else:  # Move left
            new_col = vehicle_pos[0][1] - 1
            return (new_col >= 0 and 
                   board[vehicle_pos[0][0]][new_col] == '.')
    else:  # Vertical
        if direction > 0:  # Move down
            new_row = vehicle_pos[-1][0] + 1
            return (new_row < len(board) and 
                   board[new_row][vehicle_pos[0][1]] == '.')
        else:  # Move up
            new_row = vehicle_pos[0][0] - 1
            return (new_row >= 0 and 
                   board[new_row][vehicle_pos[0][1]] == '.')

def move_vehicle(board, vehicle, positions, direction):
    new_board = [list(row) for row in board]
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        for i, pos in enumerate(positions):
            new_board[pos[0]][pos[1] + direction] = vehicle
    else:
        for i, pos in enumerate(positions):
            new_board[pos[0] + direction][pos[1]] = vehicle
            
    return new_board

def solve_puzzle():
    initial_board = [
        list('.FBBJK'),
        list('.F.IJK'),
        list('xAAI.L'),
        list('CCHI.L'),
        list('EGH...'),
        list('EGDDD.')
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        vehicles = get_vehicle_positions(current_board)
        if 'A' in vehicles:
            if vehicles['A'][-1][1] == len(current_board[0]) - 1:
                return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            # Try moving right/down
            if can_move(current_board, positions, 1):
                new_board = move_vehicle(current_board, vehicle, positions, 1)
                board_str = board_to_string(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    queue.append((new_board, moves + [f"{vehicle}+1"]))
            
            # Try moving left/up
            if can_move(current_board, positions, -1):
                new_board = move_vehicle(current_board, vehicle, positions, -1)
                board_str = board_to_string(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    queue.append((new_board, moves + [f"{vehicle}-1"]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")