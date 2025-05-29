from collections import deque

def get_red_car_exit_distance(board):
    # Find red car position and calculate distance to exit
    for i, row in enumerate(board):
        if 'A' in row:
            red_row = i
            red_right = max(j for j, c in enumerate(row) if c == 'A')
            return len(board[0]) - 1 - red_right
    return float('inf')

def get_blocking_cars(board, red_row):
    # Get cars blocking the red car's path to exit
    blocking = set()
    red_right = max(j for j, c in enumerate(board[red_row]) if c == 'A')
    for j in range(red_right + 1, len(board[0])):
        if board[red_row][j] != '.' and board[red_row][j] != 'x':
            blocking.add(board[red_row][j])
    return blocking

def get_car_orientation(board, car):
    positions = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == car]
    return 'H' if positions[0][0] == positions[1][0] else 'V'

def can_move(board, car, direction):
    positions = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == car]
    is_horizontal = positions[0][0] == positions[1][0]
    
    if is_horizontal:
        row = positions[0][0]
        if direction > 0:  # right
            rightmost = max(p[1] for p in positions)
            return rightmost + 1 < len(board[0]) and board[row][rightmost + 1] == '.'
        else:  # left
            leftmost = min(p[1] for p in positions)
            return leftmost > 0 and board[row][leftmost - 1] == '.'
    else:
        col = positions[0][1]
        if direction > 0:  # down
            bottommost = max(p[0] for p in positions)
            return bottommost + 1 < len(board) and board[bottommost + 1][col] == '.'
        else:  # up
            topmost = min(p[0] for p in positions)
            return topmost > 0 and board[topmost - 1][col] == '.'

def move_car(board, car, direction):
    new_board = [list(row) for row in board]
    positions = [(i, j) for i in range(len(board)) for j in range(len(board[i])) if board[i][j] == car]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "BBBCC.",
        "..JKDD",
        "AAJKLM",
        "I.EELM",
        "IFF..N",
        "GGHHHN"
    ]
    
    # Find red car row
    red_row = next(i for i, row in enumerate(initial_board) if 'A' in row)
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        if get_red_car_exit_distance(current_board) == 0:
            return moves
        
        # Get blocking cars first
        blocking_cars = get_blocking_cars(current_board, red_row)
        
        # Try moving blocking cars first
        for car in blocking_cars:
            for direction in [-1, 1]:
                if can_move(current_board, car, direction):
                    new_board = move_car(current_board, car, direction)
                    board_tuple = tuple(new_board)
                    
                    if board_tuple not in visited:
                        visited.add(board_tuple)
                        move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                        new_moves = moves + [move_str]
                        queue.appendleft((new_board, new_moves))  # Priority to blocking cars
        
        # Try moving other cars
        for i in range(len(current_board)):
            for j in range(len(current_board[0])):
                car = current_board[i][j]
                if car not in ('.', 'x') and car not in blocking_cars:
                    for direction in [-1, 1]:
                        if can_move(current_board, car, direction):
                            new_board = move_car(current_board, car, direction)
                            board_tuple = tuple(new_board)
                            
                            if board_tuple not in visited:
                                visited.add(board_tuple)
                                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                                new_moves = moves + [move_str]
                                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")