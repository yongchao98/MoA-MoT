from collections import deque

def get_red_car_path(board):
    # Find red car position and target path
    for i, row in enumerate(board):
        if 'A' in row:
            red_row = i
            red_col = row.index('A')
            path = set((red_row, j) for j in range(red_col + 2, len(row)))
            return red_row, red_col, path
    return None

def get_blocking_cars(board):
    red_row, red_col, path = get_red_car_path(board)
    blocking = set()
    
    # Find cars blocking the path
    for pos_row, pos_col in path:
        if board[pos_row][pos_col].isalpha():
            car = board[pos_row][pos_col]
            if car != '.' and car != 'x':
                blocking.add(car)
    
    return blocking

def get_moves(board, blocking_cars):
    moves = []
    height = len(board)
    width = len(board[0])
    
    # First priority: moves for blocking cars
    for i in range(height):
        for j in range(width):
            if board[i][j].isalpha() and board[i][j] != 'x':
                car = board[i][j]
                
                # Find if car is horizontal or vertical
                is_horizontal = (j + 1 < width and board[i][j] == board[i][j + 1])
                
                if is_horizontal:
                    # Find car length
                    car_length = 1
                    while j + car_length < width and board[i][j + car_length] == car:
                        car_length += 1
                    
                    # Try moving left
                    if j > 0 and board[i][j - 1] == '.':
                        priority = 1 if car in blocking_cars else 2
                        moves.append((priority, (car, -1)))
                    
                    # Try moving right
                    if j + car_length < width and board[i][j + car_length] == '.':
                        priority = 1 if car in blocking_cars else 2
                        moves.append((priority, (car, 1)))
                
                elif i + 1 < height and board[i + 1][j] == car:  # vertical car
                    # Find car length
                    car_length = 1
                    while i + car_length < height and board[i + car_length][j] == car:
                        car_length += 1
                    
                    # Try moving up
                    if i > 0 and board[i - 1][j] == '.':
                        priority = 1 if car in blocking_cars else 2
                        moves.append((priority, (car, -1)))
                    
                    # Try moving down
                    if i + car_length < height and board[i + car_length][j] == '.':
                        priority = 1 if car in blocking_cars else 2
                        moves.append((priority, (car, 1)))
    
    return [move[1] for move in sorted(moves)]

def apply_move(board, move):
    car, direction = move
    new_board = [list(row) for row in board]
    height = len(board)
    width = len(board[0])
    
    # Find car positions
    positions = []
    for i in range(height):
        for j in range(width):
            if board[i][j] == car:
                positions.append((i, j))
    
    positions.sort()
    
    # Determine if car is horizontal
    is_horizontal = positions[0][0] == positions[-1][0]
    
    if is_horizontal:
        row = positions[0][0]
        if direction < 0:  # Move left
            new_board[row][positions[0][1] - 1] = car
            new_board[row][positions[-1][1]] = '.'
        else:  # Move right
            new_board[row][positions[-1][1] + 1] = car
            new_board[row][positions[0][1]] = '.'
    else:  # Vertical car
        col = positions[0][1]
        if direction < 0:  # Move up
            new_board[positions[0][0] - 1][col] = car
            new_board[positions[-1][0]][col] = '.'
        else:  # Move down
            new_board[positions[-1][0] + 1][col] = car
            new_board[positions[0][0]][col] = '.'
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    seen = {tuple(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        for row in current_board:
            if 'A' in row and row.index('A') == len(row) - 2:
                return moves
        
        # Get blocking cars and possible moves
        blocking_cars = get_blocking_cars(current_board)
        possible_moves = get_moves(current_board, blocking_cars)
        
        for move in possible_moves:
            new_board = apply_move(current_board, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + [f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"]
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