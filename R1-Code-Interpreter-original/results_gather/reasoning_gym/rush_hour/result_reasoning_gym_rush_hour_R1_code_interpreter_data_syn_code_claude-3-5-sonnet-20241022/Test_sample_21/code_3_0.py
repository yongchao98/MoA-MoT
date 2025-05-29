from collections import deque

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha() and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_blocking_cars(board, cars):
    if 'A' not in cars:
        return set()
    a_row = cars['A'][0][0]
    a_right = max(pos[1] for pos in cars['A'])
    blocking = set()
    for col in range(a_right + 1, 6):
        if board[a_row][col] != '.' and board[a_row][col] != 'x':
            blocking.add(board[a_row][col])
    return blocking

def get_valid_moves(board, cars, blocking_cars):
    moves = []
    # First priority: Move blocking cars
    for car in blocking_cars:
        if car not in cars:
            continue
        positions = cars[car]
        is_horiz = is_horizontal(positions)
        
        if not is_horiz:  # Vertical cars blocking A
            min_row = min(pos[0] for pos in positions)
            max_row = max(pos[0] for pos in positions)
            # Prefer moving up if possible
            if min_row > 0 and board[min_row - 1][positions[0][1]] == '.':
                moves.append((car, -1))
            if max_row < 5 and board[max_row + 1][positions[0][1]] == '.':
                moves.append((car, 1))

    # Second priority: Move car A if path is clear
    if 'A' in cars:
        positions = cars['A']
        max_col = max(pos[1] for pos in positions)
        if max_col < 5 and board[positions[0][0]][max_col + 1] == '.':
            moves.append(('A', 1))

    # Last priority: Other cars that might help
    for car, positions in cars.items():
        if car in blocking_cars or car == 'A':
            continue
        is_horiz = is_horizontal(positions)
        
        if is_horiz:
            min_col = min(pos[1] for pos in positions)
            max_col = max(pos[1] for pos in positions)
            if min_col > 0 and board[positions[0][0]][min_col - 1] == '.':
                moves.append((car, -1))
            if max_col < 5 and board[positions[0][0]][max_col + 1] == '.':
                moves.append((car, 1))
        else:
            min_row = min(pos[0] for pos in positions)
            max_row = max(pos[0] for pos in positions)
            if min_row > 0 and board[min_row - 1][positions[0][1]] == '.':
                moves.append((car, -1))
            if max_row < 5 and board[max_row + 1][positions[0][1]] == '.':
                moves.append((car, 1))
    
    return moves

def make_move(board, car, direction):
    board = [list(row) for row in board]
    positions = []
    for i in range(6):
        for j in range(6):
            if board[i][j] == car:
                positions.append((i, j))
    
    is_horiz = is_horizontal(positions)
    for pos in positions:
        board[pos[0]][pos[1]] = '.'
    
    for pos in positions:
        new_row = pos[0] + (0 if is_horiz else direction)
        new_col = pos[1] + (direction if is_horiz else 0)
        board[new_row][new_col] = car
    
    return [''.join(row) for row in board]

def solve_puzzle():
    initial_board = [
        '...F..',
        '..xFG.',
        'AAEFG.',
        '..E.H.',
        '.DBBH.',
        '.DCCx.'
    ]
    
    queue = deque([(initial_board, [])])
    visited = {''.join(initial_board)}
    max_moves = 10  # Limit solution length
    
    while queue:
        current_board, moves = queue.popleft()
        
        if len(moves) >= max_moves:
            continue
            
        cars = get_car_positions(current_board)
        blocking_cars = get_blocking_cars(current_board, cars)
        
        # Check if solved
        if 'A' in cars and max(pos[1] for pos in cars['A']) == 5:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars, blocking_cars)
        for car, direction in valid_moves:
            new_board = make_move(current_board, car, direction)
            board_key = ''.join(new_board)
            
            if board_key not in visited:
                visited.add(board_key)
                new_moves = moves + [f"{car}{'+' if direction > 0 else ''}{direction}"]
                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")