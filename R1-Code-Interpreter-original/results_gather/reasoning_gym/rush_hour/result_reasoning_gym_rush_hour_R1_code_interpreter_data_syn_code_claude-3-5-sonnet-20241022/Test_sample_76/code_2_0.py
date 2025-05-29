from collections import deque

def get_car_positions(board):
    cars = {}
    height = len(board)
    width = len(board[0])
    for i in range(height):
        for j in range(width):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, car_positions, direction):
    height = len(board)
    width = len(board[0])
    if is_horizontal(car_positions):
        if direction > 0:  # Move right
            last_pos = max(p[1] for p in car_positions)
            return last_pos + 1 < width and board[car_positions[0][0]][last_pos + 1] == '.'
        else:  # Move left
            first_pos = min(p[1] for p in car_positions)
            return first_pos - 1 >= 0 and board[car_positions[0][0]][first_pos - 1] == '.'
    else:  # Vertical
        if direction > 0:  # Move down
            last_pos = max(p[0] for p in car_positions)
            return last_pos + 1 < height and board[last_pos + 1][car_positions[0][1]] == '.'
        else:  # Move up
            first_pos = min(p[0] for p in car_positions)
            return first_pos - 1 >= 0 and board[first_pos - 1][car_positions[0][1]] == '.'

def move_car(board, car, positions, direction):
    new_board = []
    for row in board:
        new_board.append(list(row))
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
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

def blocks_red_car(board, car_positions):
    # Check if the car blocks the path of the red car (AA)
    red_car_row = None
    for i in range(len(board)):
        if 'A' in board[i]:
            red_car_row = i
            break
    
    if red_car_row is None:
        return False
    
    return any(pos[0] == red_car_row and pos[1] > max(j for i, j in car_positions if board[i][j] == 'A')
               for pos in car_positions)

def solve_puzzle():
    initial_board = [
        "FBBBCC",
        "F.G.DD",
        "AAG.HI",
        "..G.HI",
        ".....J",
        "EEE..J"
    ]
    
    queue = deque([(initial_board, [])])
    visited = {''.join(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if red car (AA) reached the exit
        if 'A' in cars and max(p[1] for p in cars['A']) == len(current_board[0]) - 1:
            return moves
        
        # Prioritize moving blocking cars
        blocking_moves = []
        other_moves = []
        
        for car, positions in cars.items():
            if car == 'A':
                continue
                
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = ''.join(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        move = (new_board, moves + [f"{car}{'+' if direction > 0 else '-'}1"])
                        if blocks_red_car(current_board, positions):
                            blocking_moves.append(move)
                        else:
                            other_moves.append(move)
        
        # Add red car moves
        if 'A' in cars:
            for direction in [-1, 1]:
                if can_move(current_board, cars['A'], direction):
                    new_board = move_car(current_board, 'A', cars['A'], direction)
                    board_str = ''.join(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        blocking_moves.append((new_board, moves + [f"A{'+' if direction > 0 else '-'}1"]))
        
        # Add moves to queue with priority
        for move in blocking_moves:
            queue.append(move)
        for move in other_moves:
            queue.append(move)
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")