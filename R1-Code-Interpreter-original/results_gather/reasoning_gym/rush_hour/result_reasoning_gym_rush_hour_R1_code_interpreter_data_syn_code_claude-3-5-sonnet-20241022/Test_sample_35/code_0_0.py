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
        # Check if vehicle is horizontal
        if is_horizontal(positions):
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
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if is_horizontal(positions):
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = vehicle
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        '..BBB.',
        '...FG.',
        'AA.FG.',
        'E..FH.',
        'ECCCH.',
        'E.DD..'
    ]
    
    seen = set(''.join(initial_board))
    queue = deque([(initial_board, [])])
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (AA is at the exit)
        aa_positions = vehicles['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        # Try all possible moves
        for move in get_valid_moves(current_board, vehicles):
            new_board = apply_move(current_board, vehicles, move)
            board_str = ''.join(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"]
                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")