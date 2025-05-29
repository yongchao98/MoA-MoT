from collections import deque
import copy

# Initial board state
initial_board = [
    '....H.',
    '.xBBH.',
    '.AA.I.',
    'CCFGI.',
    'E.FGDD',
    'E.x...'
]

def print_board(board):
    for row in board:
        print(row)

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def make_move(board, vehicle, positions, direction):
    board = [list(row) for row in board]
    # Clear current positions
    for i, j in positions:
        board[i][j] = '.'
    # Set new positions
    if is_horizontal(positions):
        offset = 1 if direction > 0 else -1
        for i, j in positions:
            board[i][j + direction] = vehicle
    else:
        offset = 1 if direction > 0 else -1
        for i, j in positions:
            board[i + direction][j] = vehicle
    return [''.join(row) for row in board]

def can_move(board, positions, direction):
    if is_horizontal(positions):
        row = positions[0][0]
        if direction > 0:  # Try move right
            col = max(pos[1] for pos in positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # Try move left
            col = min(pos[1] for pos in positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = positions[0][1]
        if direction > 0:  # Try move down
            row = max(pos[0] for pos in positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # Try move up
            row = min(pos[0] for pos in positions) - 1
            return row >= 0 and board[row][col] == '.'

def is_solved(board):
    # Check if AA is at the exit
    for i in range(len(board)):
        if 'AA' in board[i] and board[i].index('A') + 2 == len(board[i]):
            return True
    return False

def board_to_string(board):
    return ''.join(board)

def solve_puzzle():
    start = initial_board
    queue = deque([(start, [])])
    seen = {board_to_string(start)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = make_move(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")