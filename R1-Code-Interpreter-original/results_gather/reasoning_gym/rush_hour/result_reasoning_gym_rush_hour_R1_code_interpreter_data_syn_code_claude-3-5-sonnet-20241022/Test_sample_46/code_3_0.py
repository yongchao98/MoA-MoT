from collections import deque

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j].isalpha() and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return {car: sorted(positions) for car, positions in cars.items()}

def get_moves(board):
    moves = []
    height = len(board)
    width = len(board[0])
    
    # Find all cars and their positions
    cars = get_car_positions(board)
    
    # For each car, try moving it
    for car, positions in cars.items():
        # Check if car is horizontal
        if positions[0][0] == positions[-1][0]:
            row = positions[0][0]
            # Try moving left
            if positions[0][1] > 0 and board[row][positions[0][1]-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if positions[-1][1] < width-1 and board[row][positions[-1][1]+1] == '.':
                moves.append((car, 1))
        # Car is vertical
        else:
            col = positions[0][1]
            # Try moving up
            if positions[0][0] > 0 and board[positions[0][0]-1][col] == '.':
                moves.append((car, -1))
            # Try moving down
            if positions[-1][0] < height-1 and board[positions[-1][0]+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, move):
    car, direction = move
    new_board = [list(row) for row in board]
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == car:
                positions.append((i, j))
    
    positions = sorted(positions)
    # Check if car is horizontal
    if positions[0][0] == positions[-1][0]:
        row = positions[0][0]
        if direction < 0:  # Move left
            new_board[row][positions[0][1]-1] = car
            new_board[row][positions[-1][1]] = '.'
        else:  # Move right
            new_board[row][positions[-1][1]+1] = car
            new_board[row][positions[0][1]] = '.'
    else:  # Car is vertical
        col = positions[0][1]
        if direction < 0:  # Move up
            new_board[positions[0][0]-1][col] = car
            new_board[positions[-1][0]][col] = '.'
        else:  # Move down
            new_board[positions[-1][0]+1][col] = car
            new_board[positions[0][0]][col] = '.'
    
    return [''.join(row) for row in new_board]

def is_solved(board):
    # Check if red car (AA) can exit
    for row in board:
        if 'A' in row:
            return row.index('A') == len(row)-2
    return False

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
            
        # Get all possible moves
        possible_moves = get_moves(current_board)
        
        # Try each move
        for move in possible_moves:
            new_board = apply_move(current_board, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                move_notation = f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"
                new_moves = moves + [move_notation]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
board = [
    "xBBB.I",
    ".CC..I",
    "..GAAI",
    "..GHDD",
    "..GHEE",
    ".FFx.."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")