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

def get_possible_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            # Try move left
            if left >= 0 and board[row][left] == '.':
                moves.append((vehicle, -1))
            # Try move right
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            # Try move up
            if top >= 0 and board[top][col] == '.':
                moves.append((vehicle, -1))
            # Try move down
            if bottom < len(board) and board[bottom][col] == '.':
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

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [], [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, path, moves_list = queue.popleft()
        
        # Check if solved (AA is at the exit)
        vehicles = get_vehicle_positions(current_board)
        if any(pos[1] == len(current_board[0])-1 for pos in vehicles['A']):
            return moves_list
        
        # Try all possible moves
        possible_moves = get_possible_moves(current_board, vehicles)
        for move in possible_moves:
            new_board, new_positions = apply_move(current_board, vehicles, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                vehicle, direction = move
                move_str = f"{vehicle}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, path + [new_board], moves_list + [move_str]))
    
    return None

# Initial board
board = [
    "BBH.CC",
    ".GH..K",
    ".GAAJK",
    "DD.IJL",
    "...IxL",
    ".EEFFL"
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")