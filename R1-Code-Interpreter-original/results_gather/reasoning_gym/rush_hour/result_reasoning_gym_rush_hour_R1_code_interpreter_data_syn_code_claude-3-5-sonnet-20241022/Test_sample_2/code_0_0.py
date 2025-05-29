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
            row = positions[0][0]
            left_pos = min(p[1] for p in positions) - 1
            right_pos = max(p[1] for p in positions) + 1
            
            # Try moving left
            if left_pos >= 0 and board[row][left_pos] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            if right_pos < len(board[0]) and board[row][right_pos] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top_pos = min(p[0] for p in positions) - 1
            bottom_pos = max(p[0] for p in positions) + 1
            
            # Try moving up
            if top_pos >= 0 and board[top_pos][col] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            if bottom_pos < len(board) and board[bottom_pos][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(board, vehicles, vehicle, direction):
    new_board = [list(row) for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "BBH.CC",
        "G.H.JK",
        "G.AAJK",
        "DD.IxL",
        "EE.I.L",
        "FFF..x"
    ]
    
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        vehicles = get_vehicle_positions(current_board)
        if any(pos[0] == 2 and pos[1] == 5 for pos in vehicles['A']):
            return moves
        
        # Get and apply valid moves
        valid_moves = get_valid_moves(current_board, vehicles)
        for vehicle, direction in valid_moves:
            new_board = apply_move(current_board, vehicles, vehicle, direction)
            
            if tuple(new_board) not in seen:
                seen.add(tuple(new_board))
                new_moves = moves + [(vehicle, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format output
solution = solve_puzzle()
if solution:
    formatted_moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}{abs(direction)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")