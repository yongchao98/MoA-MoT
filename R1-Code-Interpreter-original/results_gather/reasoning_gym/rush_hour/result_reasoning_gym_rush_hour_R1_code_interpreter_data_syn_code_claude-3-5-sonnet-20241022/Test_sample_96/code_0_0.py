from collections import deque
import copy

def create_board():
    board = [
        ['B', 'B', 'C', 'C', '.', 'x'],
        ['D', 'D', 'J', 'E', 'E', 'M'],
        ['.', '.', 'J', 'A', 'A', 'M'],
        ['x', '.', 'K', 'F', 'F', 'M'],
        ['.', '.', 'K', 'L', 'G', 'G'],
        ['.', 'H', 'H', 'L', 'I', 'I']
    ]
    return board

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
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
        row = vehicle_pos[0][0]
        if direction > 0:  # move right
            col = max(pos[1] for pos in vehicle_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # move left
            col = min(pos[1] for pos in vehicle_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = vehicle_pos[0][1]
        if direction > 0:  # move down
            row = max(pos[0] for pos in vehicle_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # move up
            row = min(pos[0] for pos in vehicle_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
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
    initial_board = create_board()
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(board_to_string(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car 'A' reaches exit)
        if any(pos == (2, 4) for pos in vehicles['A']):
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")