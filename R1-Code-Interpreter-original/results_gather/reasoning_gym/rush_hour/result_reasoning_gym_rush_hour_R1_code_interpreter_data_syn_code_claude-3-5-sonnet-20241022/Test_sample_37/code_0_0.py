from collections import deque
import copy

def create_initial_board():
    board = [
        ['B', 'B', 'C', 'C', 'C', '.'],
        ['.', '.', '.', 'G', 'H', '.'],
        ['A', 'A', '.', 'G', 'H', '.'],
        ['.', '.', 'D', 'D', 'H', '.'],
        ['.', 'E', 'E', 'E', '.', '.'],
        ['.', '.', 'F', 'F', 'F', '.']
    ]
    return board

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.':
                car = board[i][j]
                if car not in cars:
                    cars[car] = [(i, j)]
                else:
                    cars[car].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in car_positions) + 1
            return col < 6 and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in car_positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in car_positions) + 1
            return row < 6 and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in car_positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [row[:] for row in board]
    # Clear current positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = car
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = car
            
    return new_board

def solve_puzzle():
    initial_board = create_initial_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) is at exit
        cars = get_car_positions(current_board)
        if any(pos[0] == 2 and pos[1] == 4 for pos in cars['A']):
            return moves
        
        # Try all possible moves for each car
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")