from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[0])):
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
        if direction > 0:  # Move right
            last_pos = max(p[1] for p in car_positions)
            return last_pos + 1 < len(board[0]) and board[car_positions[0][0]][last_pos + 1] == '.'
        else:  # Move left
            first_pos = min(p[1] for p in car_positions)
            return first_pos - 1 >= 0 and board[car_positions[0][0]][first_pos - 1] == '.'
    else:  # Vertical
        if direction > 0:  # Move down
            last_pos = max(p[0] for p in car_positions)
            return last_pos + 1 < len(board) and board[last_pos + 1][car_positions[0][1]] == '.'
        else:  # Move up
            first_pos = min(p[0] for p in car_positions)
            return first_pos - 1 >= 0 and board[first_pos - 1][car_positions[0][1]] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
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

def solve_puzzle():
    initial_board = [
        "FBBBCC",
        "F.G.DD",
        "AAG.HI",
        "..G.HI",
        "....J",
        "EEE..J"
    ]
    
    queue = deque([(initial_board, [])])
    visited = set([''.join(initial_board)])
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) reached the exit
        cars = get_car_positions(current_board)
        if max(p[1] for p in cars['A']) == len(current_board[0]) - 1:
            return moves
        
        # Try moving each car
        for car, positions in cars.items():
            # Try moving left/up (-1) and right/down (+1)
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = ''.join(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")