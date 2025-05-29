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
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < 5 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < 5 and board[bottommost+1][positions[0][1]] == '.':
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
    
    return [''.join(row) for row in new_board], {k: v if k != vehicle else new_positions for k, v in vehicles.items()}

def board_to_string(board):
    return '\n'.join(board)

def solve_puzzle(initial_board):
    initial_board = initial_board.strip().split('\n')
    vehicles = get_vehicle_positions(initial_board)
    
    queue = deque([(initial_board, vehicles, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if solved (red car AA can reach exit)
        aa_positions = current_vehicles['A']
        rightmost_aa = max(p[1] for p in aa_positions)
        if rightmost_aa == 4 and all(current_board[aa_positions[0][0]][j] == '.' for j in range(rightmost_aa + 1, 6)):
            return moves
        
        # Try all possible moves
        possible_moves = get_possible_moves(current_board, current_vehicles)
        for move in possible_moves:
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            board_str = board_to_string(new_board)
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"]
                queue.append((new_board, new_vehicles, new_moves))
    
    return None

# Initial board
initial_board = """G.IBBx
G.IJCC
HAAJ.K
HDDD.K
EEFF.K
......"""

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")