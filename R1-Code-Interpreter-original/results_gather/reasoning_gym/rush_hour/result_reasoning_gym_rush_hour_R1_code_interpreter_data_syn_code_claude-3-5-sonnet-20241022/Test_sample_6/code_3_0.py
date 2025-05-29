from collections import deque
import heapq

def create_board():
    return [
        ['B', 'B', 'H', '.', 'K', '.'],
        ['.', 'G', 'H', '.', 'K', 'L'],
        ['.', 'G', 'A', 'A', 'K', 'L'],
        ['C', 'C', 'I', '.', '.', 'x'],
        ['F', '.', 'I', 'J', 'D', 'D'],
        ['F', 'E', 'E', 'J', '.', 'x']
    ]

def get_blocking_vehicles(board):
    # Get vehicles blocking the path of car A to the exit
    blocking = set()
    row = 2  # Row where car A is located
    a_pos = None
    for j in range(6):
        if board[row][j] == 'A':
            a_pos = j
            break
    
    if a_pos is not None:
        for j in range(a_pos + 2, 5):  # Check from right of A to exit
            if board[row][j] != '.' and board[row][j] != 'x':
                blocking.add(board[row][j])
    return blocking

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
    if is_horizontal(vehicle_pos):
        row = vehicle_pos[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in vehicle_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in vehicle_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = vehicle_pos[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in vehicle_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in vehicle_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
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

def solve_puzzle():
    initial_board = create_board()
    initial_blocking = get_blocking_vehicles(initial_board)
    
    # Priority queue: (num_moves + heuristic, num_moves, board, moves)
    pq = [(len(initial_blocking), 0, initial_board, [])]
    visited = {board_to_string(initial_board)}
    
    while pq:
        _, num_moves, current_board, moves = heapq.heappop(pq)
        
        if current_board[2][4] == 'A':
            return moves
            
        vehicles = get_vehicle_positions(current_board)
        blocking = get_blocking_vehicles(current_board)
        
        for vehicle, positions in vehicles.items():
            # Prioritize moving blocking vehicles and car A
            if vehicle not in blocking and vehicle != 'A':
                continue
                
            positions = sorted(positions)
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        new_blocking = get_blocking_vehicles(new_board)
                        priority = len(new_moves) + len(new_blocking)
                        heapq.heappush(pq, (priority, len(new_moves), new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")