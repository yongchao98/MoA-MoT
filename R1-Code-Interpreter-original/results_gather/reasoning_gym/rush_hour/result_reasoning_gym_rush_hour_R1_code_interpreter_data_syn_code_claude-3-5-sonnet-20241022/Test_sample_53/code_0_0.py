from collections import deque
import copy

def get_board_string(board):
    return '\n'.join(''.join(row) for row in board)

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
    new_board = [list(row) for row in board]
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = vehicle
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = vehicle
            
    return new_board

def solve_puzzle():
    initial_board = [
        list("..Ex.."),
        list("..EG.I"),
        list("AAFGHI"),
        list(".DF.H."),
        list(".DBBH."),
        list("..CCC.")
    ]
    
    queue = deque([(initial_board, [])])
    seen = {get_board_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        vehicles = get_vehicle_positions(current_board)
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles['A']):
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            # Try moving right/down (+1) and left/up (-1)
            for direction in [1, -1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_string = get_board_string(new_board)
                    
                    if board_string not in seen:
                        seen.add(board_string)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")