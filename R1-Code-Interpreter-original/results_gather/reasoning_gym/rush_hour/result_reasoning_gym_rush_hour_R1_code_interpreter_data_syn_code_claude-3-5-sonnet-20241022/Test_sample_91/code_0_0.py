from collections import deque
import copy

def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < 5 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((car, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < 5 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        new_pos = (pos[0], pos[1] + direction) if horizontal else (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    start = tuple(initial_board)
    visited = {start}
    queue = deque([(start, [])])
    
    while queue:
        current_board, path = queue.popleft()
        cars = get_car_info(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        if 'A' in cars and any(pos[1] == 4 for pos in cars['A']):
            return path
        
        # Try all possible moves
        for move in get_moves(current_board, cars):
            new_board = apply_move(current_board, cars, move)
            new_board_tuple = tuple(new_board)
            
            if new_board_tuple not in visited:
                visited.add(new_board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, path + [move_str]))
    
    return None

# Initial board
board = [
    "BBCCC.",
    "....HI",
    "..AAHI",
    "...GDD",
    "...G..",
    ".xEEFF"
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")