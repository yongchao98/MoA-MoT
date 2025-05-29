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

def get_possible_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            if left >= 0 and board[row][left] == '.':
                moves.append((vehicle, -1))
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            up = min(p[0] for p in positions) - 1
            down = max(p[0] for p in positions) + 1
            
            if up >= 0 and board[up][col] == '.':
                moves.append((vehicle, -1))
            if down < len(board) and board[down][col] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    positions = vehicles[vehicle]
    new_board = [row[:] for row in board]
    new_vehicles = copy.deepcopy(vehicles)
    
    horizontal = is_horizontal(positions)
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    new_vehicles[vehicle] = new_positions
    return new_board, new_vehicles

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = get_vehicle_positions(board)
    
    queue = deque([(board, vehicles, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if red car (AA) can exit
        red_car = current_vehicles['A']
        if max(p[1] for p in red_car) == len(current_board[0]) - 2:
            return moves
        
        possible_moves = get_possible_moves(current_board, current_vehicles)
        for move in possible_moves:
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(move[0], move[1])]
                queue.append((new_board, new_vehicles, new_moves))
    
    return None

# Initial board
board = """BBBCCK
HDDJ.K
HAAJ.L
..IEEL
FFI...
GGG..x"""

solution = solve_puzzle(board)
if solution:
    move_str = ''
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        move_str += f"{vehicle}{sign}{abs(direction)} "
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")