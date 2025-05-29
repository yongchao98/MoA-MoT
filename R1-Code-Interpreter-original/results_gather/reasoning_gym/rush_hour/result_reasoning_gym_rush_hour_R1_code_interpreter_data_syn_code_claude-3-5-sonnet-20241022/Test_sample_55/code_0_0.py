from collections import deque
import copy

def create_board():
    board = [
        ['E','B','B','.','C','C'],
        ['E','.','F','.','.','H'],
        ['A','A','F','.','G','H'],
        ['.','.','.','.','G','I'],
        ['.','x','D','D','G','I'],
        ['.','.','.','.','.','I']
    ]
    return board

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i,j)]
                else:
                    vehicles[board[i][j]].append((i,j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, vehicle_pos, direction):
    if is_horizontal(vehicle_pos):
        if direction > 0:  # moving right
            new_col = vehicle_pos[-1][1] + 1
            return new_col < 6 and board[vehicle_pos[0][0]][new_col] == '.'
        else:  # moving left
            new_col = vehicle_pos[0][1] - 1
            return new_col >= 0 and board[vehicle_pos[0][0]][new_col] == '.'
    else:  # vertical
        if direction > 0:  # moving down
            new_row = vehicle_pos[-1][0] + 1
            return new_row < 6 and board[new_row][vehicle_pos[0][1]] == '.'
        else:  # moving up
            new_row = vehicle_pos[0][0] - 1
            return new_row >= 0 and board[new_row][vehicle_pos[0][1]] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = copy.deepcopy(board)
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    new_positions = []
    for pos in positions:
        if is_horizontal(positions):
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    return new_board, new_positions

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [], {})])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, vehicle_positions = queue.popleft()
        
        if not vehicle_positions:
            vehicle_positions = get_vehicle_positions(current_board)
            
        # Check if solved (red car AA at exit)
        if 'A' in vehicle_positions and vehicle_positions['A'][-1][1] == 4:
            return moves
            
        # Try moving each vehicle
        for vehicle, positions in vehicle_positions.items():
            # Try moving in both directions
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board, new_positions = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_vehicle_positions = copy.deepcopy(vehicle_positions)
                        new_vehicle_positions[vehicle] = new_positions
                        move_str = f"{vehicle}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str], new_vehicle_positions))

    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")