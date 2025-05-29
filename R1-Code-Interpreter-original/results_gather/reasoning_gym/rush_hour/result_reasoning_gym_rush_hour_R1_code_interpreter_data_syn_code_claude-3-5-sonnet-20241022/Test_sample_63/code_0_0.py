from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def get_orientation(positions):
    if positions[0][0] == positions[1][0]:  # same row
        return 'horizontal'
    return 'vertical'

def is_valid_move(board, car_positions, direction, distance):
    orientation = get_orientation(car_positions)
    
    if orientation == 'horizontal':
        row = car_positions[0][0]
        if direction > 0:  # moving right
            rightmost = max(pos[1] for pos in car_positions)
            return (rightmost + direction < 6 and 
                    all(board[row][rightmost + i + 1] == '.' for i in range(distance)))
        else:  # moving left
            leftmost = min(pos[1] for pos in car_positions)
            return (leftmost + direction >= 0 and 
                    all(board[row][leftmost + i] == '.' for i in range(direction, 0)))
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # moving down
            bottommost = max(pos[0] for pos in car_positions)
            return (bottommost + direction < 6 and 
                    all(board[bottommost + i + 1][col] == '.' for i in range(distance)))
        else:  # moving up
            topmost = min(pos[0] for pos in car_positions)
            return (topmost + direction >= 0 and 
                    all(board[topmost + i][col] == '.' for i in range(direction, 0)))

def make_move(board, car, car_positions, direction):
    new_board = [list(row) for row in board]
    orientation = get_orientation(car_positions)
    
    # Clear current positions
    for pos in car_positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    new_positions = []
    for pos in car_positions:
        if orientation == 'horizontal':
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (red car AA at exit)
        if 'A' in current_board[2] and current_board[2].rindex('A') == 5:
            return moves
        
        cars = get_car_positions(current_board)
        
        for car, positions in cars.items():
            orientation = get_orientation(positions)
            
            # Try moves in both directions
            for direction in [-1, 1]:
                if is_valid_move(current_board, positions, direction, 1):
                    new_board = make_move(current_board, car, positions, direction)
                    move = f"{car}{'+' if direction > 0 else '-'}1"
                    
                    if tuple(new_board) not in seen:
                        seen.add(tuple(new_board))
                        queue.append((new_board, moves + [move]))
    
    return None

# Initial board
initial_board = [
    '.BBBKM',
    'CC.IKM',
    'AA.ILM',
    'GDDJL.',
    'G.HJEE',
    'FFH...'
]

solution = solve_puzzle(initial_board)
if solution:
    print(' '.join(solution))