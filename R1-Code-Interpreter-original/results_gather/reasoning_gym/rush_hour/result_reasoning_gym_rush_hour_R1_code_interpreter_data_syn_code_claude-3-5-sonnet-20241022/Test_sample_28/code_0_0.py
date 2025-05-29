from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

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

def apply_move(board, vehicles, vehicle, direction):
    new_board = [list(row) for row in board]
    positions = vehicles[vehicle]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        new_pos = (pos[0], pos[1] + direction) if horizontal else (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = vehicle
        new_positions.append(new_pos)
    
    return new_board, new_positions

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [], [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, path, moves = queue.popleft()
        
        # Check if solved (red car AA at exit)
        vehicles = get_vehicle_positions(current_board)
        if any(pos[0] == 2 and pos[1] == 5 for pos in vehicles['A']):
            return moves
        
        # Try all possible moves
        for vehicle, direction in get_valid_moves(current_board, vehicles):
            new_board, new_positions = apply_move(current_board, vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [f"{vehicle}{'+' if direction > 0 else ''}{direction}"]
                queue.append((new_board, path + [(vehicle, direction)], new_moves))
    
    return None

# Initial board
initial_board = [
    list("..xBBJ"),
    list("CCCDDJ"),
    list("..AA.K"),
    list("G.HEEK"),
    list("G.HIFF"),
    list("...I..")
]

solution = solve_puzzle(initial_board)
if solution:
    print(' '.join(solution))