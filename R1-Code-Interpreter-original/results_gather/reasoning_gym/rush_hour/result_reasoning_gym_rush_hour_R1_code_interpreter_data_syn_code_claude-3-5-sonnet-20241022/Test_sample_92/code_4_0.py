def create_board():
    return [
        ['.', '.', 'x', 'B', 'B', 'B'],
        ['G', 'H', 'I', '.', 'C', 'C'],
        ['G', 'H', 'I', 'A', 'A', 'K'],
        ['D', 'D', 'D', 'J', '.', 'K'],
        ['.', '.', '.', 'J', 'E', 'E'],
        ['.', 'F', 'F', 'F', '.', 'x']
    ]

def find_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction > 0:  # right
            next_col = max(pos[1] for pos in car_positions) + 1
            return next_col < 6 and board[row][next_col] == '.'
        else:  # left
            next_col = min(pos[1] for pos in car_positions) - 1
            return next_col >= 0 and board[row][next_col] == '.'
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # down
            next_row = max(pos[0] for pos in car_positions) + 1
            return next_row < 6 and board[next_row][col] == '.'
        else:  # up
            next_row = min(pos[0] for pos in car_positions) - 1
            return next_row >= 0 and board[next_row][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [row[:] for row in board]
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        for i, j in positions:
            new_board[row][j + direction] = car
    else:
        col = positions[0][1]
        for i, j in positions:
            new_board[i + direction][col] = car
            
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_rush_hour():
    initial_board = create_board()
    visited = {board_to_string(initial_board)}
    queue = [(initial_board, [])]
    
    while queue:
        current_board, moves = queue.pop(0)
        cars = find_car_positions(current_board)
        
        # Check if red car can reach exit
        red_car = cars['A']
        red_rightmost = max(j for _, j in red_car)
        can_exit = True
        for j in range(red_rightmost + 1, 6):
            if current_board[2][j] not in '.x':
                can_exit = False
                break
        
        if can_exit:
            moves.append(f"A+{5 - red_rightmost}")
            return moves
        
        # Try all possible moves
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}1"]
                        queue.append((new_board, new_moves))
    
    return None  # No solution found

solution = solve_rush_hour()
if solution is None:
    print("<<<No solution found>>>")
else:
    print('<<<' + ' '.join(solution) + '>>>')