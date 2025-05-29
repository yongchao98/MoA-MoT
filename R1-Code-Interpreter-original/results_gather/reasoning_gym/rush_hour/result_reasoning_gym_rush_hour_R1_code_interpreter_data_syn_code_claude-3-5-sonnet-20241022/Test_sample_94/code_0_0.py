from collections import deque

def create_board(cars):
    board = [['.'] * 6 for _ in range(6)]
    board[4][3] = 'x'  # wall
    for car, (r, c, length, is_horizontal) in cars.items():
        for i in range(length):
            if is_horizontal:
                board[r][c + i] = car
            else:
                board[r + i][c] = car
    return board

def get_moves(board, cars):
    moves = []
    for car, (r, c, length, is_horizontal) in cars.items():
        if is_horizontal:
            # Try moving left
            if c > 0 and board[r][c-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if c + length < 6 and board[r][c+length] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            if r > 0 and board[r-1][c] == '.':
                moves.append((car, -1))
            # Try moving down
            if r + length < 6 and board[r+length][c] == '.':
                moves.append((car, 1))
    return moves

def apply_move(cars, car, delta):
    new_cars = dict(cars)
    r, c, length, is_horizontal = new_cars[car]
    if is_horizontal:
        new_cars[car] = (r, c + delta, length, is_horizontal)
    else:
        new_cars[car] = (r + delta, c, length, is_horizontal)
    return new_cars

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    # Initial state
    cars = {
        'A': (2, 1, 2, True),   # red car
        'B': (1, 2, 3, True),
        'C': (4, 1, 2, True),
        'D': (0, 2, 2, False),
        'E': (2, 3, 2, False),
        'F': (2, 5, 3, False)
    }
    
    start_board = create_board(cars)
    visited = {board_to_string(start_board)}
    queue = deque([(cars, [])])
    
    while queue:
        current_cars, moves = queue.popleft()
        board = create_board(current_cars)
        
        # Check if red car (A) can exit
        r, c, length, _ = current_cars['A']
        if c + length == 5:
            return moves
        
        for car, delta in get_moves(board, current_cars):
            new_cars = apply_move(current_cars, car, delta)
            new_board = create_board(new_cars)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [f"{car}{'+' if delta > 0 else ''}{delta}"]
                queue.append((new_cars, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")