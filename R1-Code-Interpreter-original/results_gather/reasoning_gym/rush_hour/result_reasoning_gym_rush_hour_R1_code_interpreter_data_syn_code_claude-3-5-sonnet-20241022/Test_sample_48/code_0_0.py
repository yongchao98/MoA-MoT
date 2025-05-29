from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        # Horizontal car
        if is_horizontal(positions):
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        # Vertical car
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((car, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((car, 1))

    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    if is_horizontal(positions):
        new_positions = [(p[0], p[1] + direction) for p in positions]
    else:
        new_positions = [(p[0] + direction, p[1]) for p in positions]
    
    # Place car in new positions
    for pos in new_positions:
        new_board[pos[0]][pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (AA reaches right edge)
        aa_positions = cars['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return moves
        
        # Try all possible moves
        for move in get_valid_moves(current_board, cars):
            new_board = apply_move(current_board, cars, move)
            new_board_tuple = tuple(new_board)
            
            if new_board_tuple not in seen:
                seen.add(new_board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = [
    "H.BBCC",
    "HDD..L",
    "IAA..L",
    "IEEE.L",
    "FFJKGG",
    "..JK.."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")