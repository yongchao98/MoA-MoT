from collections import deque

def get_red_car_position(board):
    for i, row in enumerate(board):
        if 'A' in row:
            start = row.index('A')
            end = row.rindex('A')
            return (i, start, end)
    return None

def is_blocking(car_pos, red_row, red_end):
    if isinstance(car_pos[0], int):  # Single position
        return car_pos[0] == red_row and car_pos[1] > red_end
    else:  # Multiple positions
        return any(pos[0] == red_row and pos[1] > red_end for pos in car_pos)

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return {k: sorted(v) for k, v in cars.items()}

def is_horizontal(positions):
    return all(pos[0] == positions[0][0] for pos in positions)

def get_valid_moves(board, cars):
    moves = []
    red_pos = get_red_car_position(board)
    red_row, red_start, red_end = red_pos
    
    # First priority: Try to move red car if path is clear
    if 'A' in cars:
        if red_end < 5 and board[red_row][red_end + 1] == '.':
            moves.append(('A', 1))
        if red_start > 0 and board[red_row][red_start - 1] == '.':
            moves.append(('A', -1))

    # Second priority: Move blocking cars
    blocking_cars = []
    for car, positions in cars.items():
        if car != 'A' and is_blocking(positions, red_row, red_end):
            blocking_cars.append(car)
            
    for car in blocking_cars:
        positions = cars[car]
        if is_horizontal(positions):
            # Try horizontal moves
            if positions[0][1] > 0 and board[positions[0][0]][positions[0][1] - 1] == '.':
                moves.append((car, -1))
            if positions[-1][1] < 5 and board[positions[0][0]][positions[-1][1] + 1] == '.':
                moves.append((car, 1))
        else:
            # Try vertical moves
            if positions[0][0] > 0 and board[positions[0][0] - 1][positions[0][1]] == '.':
                moves.append((car, -1))
            if positions[-1][0] < 5 and board[positions[-1][0] + 1][positions[0][1]] == '.':
                moves.append((car, 1))

    # Last priority: Move other cars
    for car, positions in cars.items():
        if car not in blocking_cars and car != 'A':
            if is_horizontal(positions):
                if positions[0][1] > 0 and board[positions[0][0]][positions[0][1] - 1] == '.':
                    moves.append((car, -1))
                if positions[-1][1] < 5 and board[positions[0][0]][positions[-1][1] + 1] == '.':
                    moves.append((car, 1))
            else:
                if positions[0][0] > 0 and board[positions[0][0] - 1][positions[0][1]] == '.':
                    moves.append((car, -1))
                if positions[-1][0] < 5 and board[positions[-1][0] + 1][positions[0][1]] == '.':
                    moves.append((car, 1))
    
    return moves

def apply_move(board, car, positions, direction):
    new_board = [list(row) for row in board]
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    is_hor = is_horizontal(positions)
    for pos in positions:
        new_row = pos[0] if is_hor else pos[0] + direction
        new_col = pos[1] + direction if is_hor else pos[1]
        new_board[new_row][new_col] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "..HBBB",
        "..HICC",
        ".AAIJK",
        ".GDDJK",
        ".GEEJL",
        "FF...L"
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    
    while queue:
        current_board, moves = queue.popleft()
        board_str = '\n'.join(current_board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        # Check if solved
        red_pos = get_red_car_position(current_board)
        if red_pos and red_pos[2] == 5:
            return moves
        
        cars = get_car_positions(current_board)
        valid_moves = get_valid_moves(current_board, cars)
        
        for car, direction in valid_moves:
            new_board = apply_move(current_board, car, cars[car], direction)
            move_str = f"{car}{'+' if direction > 0 else '-'}1"
            
            # Prioritize promising moves
            if car == 'A' and direction > 0:
                queue.appendleft((new_board, moves + [move_str]))
            else:
                queue.append((new_board, moves + [move_str]))

    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")