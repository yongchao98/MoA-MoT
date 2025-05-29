from collections import deque
import copy

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

def get_valid_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((vehicle, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    positions = vehicles[vehicle]
    new_board = [list(row) for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board], new_positions

def board_to_string(board):
    return '|'.join(board)

def solve_puzzle():
    initial_board = [
        'x.BBCC',
        'DDEELM',
        'IJAALM',
        'IJFFLN',
        'GGK..N',
        '..K.HH'
    ]
    
    queue = deque([(initial_board, [], {})])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, path, _ = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (AA is at the exit)
        aa_positions = vehicles['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return path
        
        # Try all possible moves
        for move in get_valid_moves(current_board, vehicles):
            new_board, new_positions = apply_move(current_board, vehicles, move)
            board_string = board_to_string(new_board)
            
            if board_string not in seen:
                seen.add(board_string)
                vehicle, direction = move
                move_str = f"{vehicle}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, path + [move_str], {}))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")