from collections import deque

def parse_board(board):
    return [list(row) for row in board.strip().split('\n')]

def get_blocking_vehicles(board):
    # Find red car row and position
    red_row = None
    red_right = 0
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                red_row = i
                red_right = max(red_right, j)
    
    # Get vehicles blocking the path
    blocking = set()
    for j in range(red_right + 1, len(board[0])-1):
        if board[red_row][j] not in '.x':
            blocking.add(board[red_row][j])
    return blocking

def get_vehicle_positions(board, vehicle):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == vehicle:
                positions.append((i, j))
    return positions

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, vehicle):
    positions = get_vehicle_positions(board, vehicle)
    horizontal = is_horizontal(positions)
    
    if horizontal:
        row = positions[0][0]
        min_col = min(p[1] for p in positions)
        max_col = max(p[1] for p in positions)
        return (min_col > 0 and board[row][min_col-1] == '.', 
                max_col < len(board[0])-1 and board[row][max_col+1] == '.')
    else:
        col = positions[0][1]
        min_row = min(p[0] for p in positions)
        max_row = max(p[0] for p in positions)
        return (min_row > 0 and board[min_row-1][col] == '.',
                max_row < len(board)-1 and board[max_row+1][col] == '.')

def move_vehicle(board, vehicle, direction):
    new_board = [row[:] for row in board]
    positions = get_vehicle_positions(board, vehicle)
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    for i, j in positions:
        new_i = i + (0 if horizontal else direction)
        new_j = j + (direction if horizontal else 0)
        new_board[new_i][new_j] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_rush_hour(board):
    initial_board = parse_board(board)
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car can reach exit
        red_positions = get_vehicle_positions(current_board, 'A')
        red_right = max(p[1] for p in red_positions)
        red_row = red_positions[0][0]
        can_exit = True
        
        for j in range(red_right + 1, len(current_board[0])-1):
            if current_board[red_row][j] != '.':
                can_exit = False
                break
        
        if can_exit:
            return moves + [('A', 4)]  # Move red car to exit
        
        # Get blocking vehicles and try to move them
        blocking = get_blocking_vehicles(current_board)
        for vehicle in blocking:
            left_possible, right_possible = can_move(current_board, vehicle)
            
            if left_possible:
                new_board = move_vehicle(current_board, vehicle, -1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, moves + [(vehicle, -1)]))
            
            if right_possible:
                new_board = move_vehicle(current_board, vehicle, 1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, moves + [(vehicle, 1)]))
        
        # Try moving other vehicles that might help
        for i in range(len(current_board)):
            for j in range(len(current_board[0])):
                vehicle = current_board[i][j]
                if vehicle not in '.x' and vehicle not in blocking and vehicle != 'A':
                    left_possible, right_possible = can_move(current_board, vehicle)
                    
                    if left_possible:
                        new_board = move_vehicle(current_board, vehicle, -1)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            queue.append((new_board, moves + [(vehicle, -1)]))
                    
                    if right_possible:
                        new_board = move_vehicle(current_board, vehicle, 1)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            queue.append((new_board, moves + [(vehicle, 1)]))
    
    return None

board = """HBBCCx
H.IDDD
AAIJ..
.EEJFF
.GGG.K
.....K"""

solution = solve_rush_hour(board)
if solution:
    moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        moves.append(f"{vehicle}{sign}1")
    print("<<<" + " ".join(moves) + ">>>")
else:
    print("<<<No solution found>>>")