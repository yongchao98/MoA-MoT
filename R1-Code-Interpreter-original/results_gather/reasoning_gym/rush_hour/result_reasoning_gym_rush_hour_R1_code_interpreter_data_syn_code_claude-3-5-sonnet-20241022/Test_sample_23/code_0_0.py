from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

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
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            
            # Try moving right
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            
            # Try moving down
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
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
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [], [])])
    seen = {tuple(initial_board)}
    
    while queue:
        board, path, moves = queue.popleft()
        cars = get_car_positions(board)
        
        # Check if solved (car AA is at the exit)
        if 'A' in cars and any(pos[1] == len(board[0])-1 for pos in cars['A']):
            return moves
        
        # Try all possible moves
        for move in get_valid_moves(board, cars):
            new_board, new_positions = apply_move(board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, path + [new_board], moves + [move_str]))
    
    return None

# Initial board
initial_board = [
    "..BBBI",
    "x.CCCI",
    "AA.GHI",
    "...GH.",
    "..DDx.",
    "..EEFF"
]

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")