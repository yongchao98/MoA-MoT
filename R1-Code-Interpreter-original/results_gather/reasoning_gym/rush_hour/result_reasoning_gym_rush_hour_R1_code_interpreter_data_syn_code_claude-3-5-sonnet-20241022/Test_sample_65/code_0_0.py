from collections import deque
import copy

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

def solve_puzzle(initial_board):
    target_row = 2  # Row where AA is located
    target_col = 5  # Exit column
    
    visited = set()
    queue = deque([(initial_board, [], [])])
    
    while queue:
        current_board, path, moves_list = queue.popleft()
        board_str = '\n'.join(current_board)
        
        if board_str in visited:
            continue
        
        visited.add(board_str)
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved
        if 'A' in vehicles and any(pos[0] == target_row and pos[1] == target_col-1 for pos in vehicles['A']):
            return moves_list + [('A', 1)]
        
        possible_moves = get_possible_moves(current_board, vehicles)
        
        for move in possible_moves:
            new_board, new_positions = apply_move(current_board, vehicles, move)
            vehicles[move[0]] = new_positions
            queue.append((new_board, path + [new_board], moves_list + [move]))
    
    return None

# Initial board
board = [
    'E..H..',
    'E..Hxx',
    'EAAIJ.',
    '.FGIJK',
    '.FGBBK',
    'CCDDD.'
]

solution = solve_puzzle(board)
if solution:
    # Format the solution
    formatted_solution = ''
    for vehicle, spaces in solution:
        sign = '+' if spaces > 0 else ''
        formatted_solution += f'{vehicle}{sign}{spaces} '
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("No solution found")