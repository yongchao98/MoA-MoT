from heapq import heappush, heappop

def parse_board(board):
    return [list(row) for row in board.strip().split('\n')]

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board, vehicles):
    # Calculate blocking vehicles between red car and exit
    red_car = None
    blocking = 0
    red_row = None
    red_right = 0
    
    # Find red car position
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                red_row = i
                red_right = max(red_right, j)
    
    # Count blocking vehicles
    for j in range(red_right + 1, len(board[0])-1):
        if board[red_row][j] not in '.x':
            blocking += 1
    
    return blocking * 2  # Weight for blocking vehicles

def get_vehicle_data(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            char = board[i][j]
            if char not in '.x':
                if char not in vehicles:
                    vehicles[char] = {'pos': set(), 'horizontal': None}
                vehicles[char]['pos'].add((i, j))
    
    for v in vehicles:
        pos_list = list(vehicles[v]['pos'])
        vehicles[v]['horizontal'] = pos_list[0][0] == pos_list[1][0]
    return vehicles

def get_moves(board, vehicles):
    moves = []
    for vehicle, data in vehicles.items():
        pos = list(data['pos'])
        if data['horizontal']:
            row = pos[0][0]
            min_col = min(p[1] for p in pos)
            max_col = max(p[1] for p in pos)
            
            if min_col > 0 and board[row][min_col-1] == '.':
                moves.append((vehicle, -1))
            if max_col < len(board[0])-1 and board[row][max_col+1] == '.':
                moves.append((vehicle, 1))
        else:
            col = pos[0][1]
            min_row = min(p[0] for p in pos)
            max_row = max(p[0] for p in pos)
            
            if min_row > 0 and board[min_row-1][col] == '.':
                moves.append((vehicle, -1))
            if max_row < len(board)-1 and board[max_row+1][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    new_board = [row[:] for row in board]
    pos = list(vehicles[vehicle]['pos'])
    horizontal = vehicles[vehicle]['horizontal']
    
    # Clear current positions
    for i, j in pos:
        new_board[i][j] = '.'
    
    # Add new positions
    for i, j in pos:
        new_i = i + (0 if horizontal else direction)
        new_j = j + (direction if horizontal else 0)
        new_board[new_i][new_j] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve(initial_board):
    board = parse_board(initial_board)
    initial_vehicles = get_vehicle_data(board)
    
    # Priority queue: (priority, moves_count, board_string, path)
    queue = [(0, 0, board_to_string(board), [])]
    seen = {board_to_string(board)}
    
    while queue:
        _, moves_count, board_str, path = heappop(queue)
        current_board = parse_board(board_str)
        current_vehicles = get_vehicle_data(current_board)
        
        # Check if solved
        red_pos = list(current_vehicles['A']['pos'])
        if max(p[1] for p in red_pos) == len(current_board[0])-2:
            return path
        
        for move in get_moves(current_board, current_vehicles):
            new_board = apply_move(current_board, current_vehicles, move)
            new_board_str = board_to_string(new_board)
            
            if new_board_str not in seen:
                seen.add(new_board_str)
                new_path = path + [move]
                priority = moves_count + 1 + heuristic(new_board, get_vehicle_data(new_board))
                heappush(queue, (priority, moves_count + 1, new_board_str, new_path))
    
    return None

board = """HBBCCx
H.IDDD
AAIJ..
.EEJFF
.GGG.K
.....K"""

solution = solve(board)
if solution:
    moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        moves.append(f"{vehicle}{sign}1")
    print("<<<" + " ".join(moves) + ">>>")
else:
    print("<<<No solution found>>>")