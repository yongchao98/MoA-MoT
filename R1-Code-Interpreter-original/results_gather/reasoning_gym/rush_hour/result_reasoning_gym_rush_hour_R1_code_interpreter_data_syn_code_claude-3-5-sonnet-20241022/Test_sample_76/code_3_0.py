from collections import deque

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

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction > 0:  # right
            max_col = max(p[1] for p in car_positions)
            return max_col + 1 < len(board[0]) and board[row][max_col + 1] == '.'
        else:  # left
            min_col = min(p[1] for p in car_positions)
            return min_col - 1 >= 0 and board[row][min_col - 1] == '.'
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # down
            max_row = max(p[0] for p in car_positions)
            return max_row + 1 < len(board) and board[max_row + 1][col] == '.'
        else:  # up
            min_row = min(p[0] for p in car_positions)
            return min_row - 1 >= 0 and board[min_row - 1][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    
    # Remove car from current position
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add car to new position
    if is_horizontal(positions):
        row = positions[0][0]
        min_col = min(p[1] for p in positions)
        for i in range(len(positions)):
            new_board[row][min_col + i + direction] = car
    else:
        col = positions[0][1]
        min_row = min(p[0] for p in positions)
        for i in range(len(positions)):
            new_board[min_row + i + direction][col] = car
            
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "FBBBCC",
        "F.G.DD",
        "AAG.HI",
        "..G.HI",
        ".....J",
        "EEE..J"
    ]
    
    # Priority order of cars to move (cars blocking red car's path first)
    priority_cars = ['G', 'H', 'I', 'J', 'D', 'C', 'A']
    
    queue = deque([(initial_board, [])])
    visited = {''.join(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if red car reached exit
        if 'A' in cars:
            red_car_pos = cars['A']
            if max(p[1] for p in red_car_pos) == len(current_board[0]) - 1:
                return moves
        
        # Try moves in priority order
        for car in priority_cars:
            if car not in cars:
                continue
                
            positions = cars[car]
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = ''.join(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}1"]
                        queue.append((new_board, new_moves))

    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")