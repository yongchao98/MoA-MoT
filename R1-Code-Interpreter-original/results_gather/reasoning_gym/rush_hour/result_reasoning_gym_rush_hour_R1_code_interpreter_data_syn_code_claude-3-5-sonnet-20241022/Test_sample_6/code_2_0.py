from collections import deque
import copy

def create_board():
    return [
        ['B', 'B', 'H', '.', 'K', '.'],
        ['.', 'G', 'H', '.', 'K', 'L'],
        ['.', 'G', 'A', 'A', 'K', 'L'],
        ['C', 'C', 'I', '.', '.', 'x'],
        ['F', '.', 'I', 'J', 'D', 'D'],
        ['F', 'E', 'E', 'J', '.', 'x']
    ]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return len(set(pos[0] for pos in positions)) == 1

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, vehicle_pos, direction):
    if len(vehicle_pos) < 2:
        return False
        
    if is_horizontal(vehicle_pos):
        row = vehicle_pos[0][0]
        if direction > 0:  # moving right
            col = max(pos[1] for pos in vehicle_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # moving left
            col = min(pos[1] for pos in vehicle_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = vehicle_pos[0][1]
        if direction > 0:  # moving down
            row = max(pos[0] for pos in vehicle_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # moving up
            row = min(pos[0] for pos in vehicle_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = copy.deepcopy(board)
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    if is_horizontal(positions):
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = vehicle
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = vehicle
            
    return new_board

def manhattan_distance_to_exit(positions, target_vehicle='A'):
    if target_vehicle == 'A':
        # Distance from rightmost position of A to exit (position 4)
        return abs(4 - max(pos[1] for pos in positions))
    return 0

def solve_puzzle():
    initial_board = create_board()
    queue = [(0, initial_board, [])]  # (priority, board, moves)
    visited = {board_to_string(initial_board)}
    
    while queue:
        queue.sort()  # Sort by priority
        _, current_board, moves = queue.pop(0)
        
        if current_board[2][4] == 'A':
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        
        for vehicle, positions in vehicles.items():
            positions = sorted(positions)
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        # Calculate priority based on Manhattan distance
                        priority = len(new_moves) + manhattan_distance_to_exit(get_vehicle_positions(new_board)['A'])
                        queue.append((priority, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")