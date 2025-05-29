from collections import deque

def parse_board(board_str):
    return [list(line) for line in board_str.splitlines()]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def can_move(board, positions, di, dj):
    for i, j in positions:
        ni, nj = i + di, j + dj
        if not (0 <= ni < len(board) and 0 <= nj < len(board[0])):
            return False
        if board[ni][nj] not in ['.', board[i][j]]:
            return False
    return True

def move_vehicle(board, positions, di, dj):
    new_board = [row[:] for row in board]
    vehicle = board[positions[0][0]][positions[0][1]]
    
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    new_positions = []
    for i, j in positions:
        ni, nj = i + di, j + dj
        new_board[ni][nj] = vehicle
        new_positions.append((ni, nj))
    
    return new_board, new_positions

def solve_rush_hour(board_str):
    board = parse_board(board_str)
    initial_vehicles = get_vehicle_positions(board)
    
    # Determine vehicle orientations
    orientations = {}
    for v, pos in initial_vehicles.items():
        orientations[v] = 'horizontal' if pos[0][0] == pos[1][0] else 'vertical'
    
    queue = deque([(board, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, moves = queue.popleft()
        current_vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (red car at exit)
        if max(j for i, j in current_vehicles['A']) == len(current_board[0]) - 1:
            return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in current_vehicles.items():
            if orientations[vehicle] == 'horizontal':
                # Try moving left and right
                for dj in [-1, 1]:
                    if can_move(current_board, positions, 0, dj):
                        new_board, new_positions = move_vehicle(current_board, positions, 0, dj)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            new_moves = moves + [f"{vehicle}{'+' if dj > 0 else ''}{dj}"]
                            queue.append((new_board, new_moves))
            else:
                # Try moving up and down
                for di in [-1, 1]:
                    if can_move(current_board, positions, di, 0):
                        new_board, new_positions = move_vehicle(current_board, positions, di, 0)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            new_moves = moves + [f"{vehicle}{'+' if di > 0 else ''}{di}"]
                            queue.append((new_board, new_moves))
    
    return None

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

# Initial board
board = """BBBCCM
..JDDM
AAJK..
IEEKFF
IGGKLx
HH..L."""

solution = solve_rush_hour(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")