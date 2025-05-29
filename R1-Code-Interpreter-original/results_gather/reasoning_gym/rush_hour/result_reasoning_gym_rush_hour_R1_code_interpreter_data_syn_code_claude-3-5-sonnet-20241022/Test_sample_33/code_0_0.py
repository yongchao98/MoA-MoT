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
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board], {**cars, car: new_positions}

def solve_puzzle(initial_board):
    board = initial_board.split('\n')
    cars = get_car_positions(board)
    visited = set([''.join(board)])
    queue = deque([(board, cars, [])])
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        aa_positions = current_cars['A']
        if max(p[1] for p in aa_positions) == 5:
            return moves
        
        # Try all possible moves
        valid_moves = get_valid_moves(current_board, current_cars)
        for move in valid_moves:
            new_board, new_cars = apply_move(current_board, current_cars, move)
            board_str = ''.join(new_board)
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [(move[0], move[1])]
                queue.append((new_board, new_cars, new_moves))
    
    return None

# Initial board
initial_board = """..xBBK
GCCJ.K
GAAJ..
DDIJ..
.HIEE.
.HFFx."""

solution = solve_puzzle(initial_board)
if solution:
    move_str = ''
    for car, spaces in solution:
        move_str += f'{car}{"+"+str(spaces) if spaces > 0 else str(spaces)} '
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")