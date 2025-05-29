from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n')]

def get_vehicle_info(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = {'pos': [], 'orientation': None}
                vehicles[board[i][j]]['pos'].append((i, j))
    
    for v in vehicles:
        if vehicles[v]['pos'][0][0] == vehicles[v]['pos'][1][0]:
            vehicles[v]['orientation'] = 'horizontal'
        else:
            vehicles[v]['orientation'] = 'vertical'
    return vehicles

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def make_move(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    # Set new positions
    if direction > 0:  # right/down
        for pos in positions:
            if vehicle['orientation'] == 'horizontal':
                new_board[pos[0]][pos[1] + direction] = board[pos[0]][pos[1]]
            else:
                new_board[pos[0] + direction][pos[1]] = board[pos[0]][pos[1]]
    else:  # left/up
        for pos in positions:
            if vehicle['orientation'] == 'horizontal':
                new_board[pos[0]][pos[1] + direction] = board[pos[0]][pos[1]]
            else:
                new_board[pos[0] + direction][pos[1]] = board[pos[0]][pos[1]]
    return new_board

def get_possible_moves(board, vehicles):
    moves = []
    for v_name, v_info in vehicles.items():
        if v_info['orientation'] == 'horizontal':
            # Try moving left
            leftmost = min(pos[1] for pos in v_info['pos'])
            if leftmost > 0 and board[v_info['pos'][0][0]][leftmost-1] == '.':
                moves.append((v_name, -1))
            # Try moving right
            rightmost = max(pos[1] for pos in v_info['pos'])
            if rightmost < 5 and board[v_info['pos'][0][0]][rightmost+1] == '.':
                moves.append((v_name, 1))
        else:
            # Try moving up
            topmost = min(pos[0] for pos in v_info['pos'])
            if topmost > 0 and board[topmost-1][v_info['pos'][0][1]] == '.':
                moves.append((v_name, -1))
            # Try moving down
            bottommost = max(pos[0] for pos in v_info['pos'])
            if bottommost < 5 and board[bottommost+1][v_info['pos'][0][1]] == '.':
                moves.append((v_name, 1))
    return moves

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = get_vehicle_info(board)
    
    queue = deque([(board, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA reaches right edge)
        for i in range(len(current_board)):
            if current_board[i][-1] == 'A':
                return moves
        
        # Get current vehicle positions
        vehicles = get_vehicle_info(current_board)
        possible_moves = get_possible_moves(current_board, vehicles)
        
        for v_name, direction in possible_moves:
            new_board = make_move(current_board, vehicles[v_name], vehicles[v_name]['pos'], direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(v_name, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
board = """BBJCC.
DDJ.L.
HIAAL.
HIEELx
H..KFF
GG.K.."""

solution = solve_puzzle(board)
if solution:
    move_str = ""
    for move in solution:
        direction = "+" if move[1] > 0 else "-"
        move_str += f"{move[0]}{direction}{abs(move[1])} "
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")