from collections import deque
import copy

# Define the initial board state
initial_board = [
    ['.', '.', 'B', 'B', 'B', 'I'],
    ['.', '.', '.', '.', 'H', 'I'],
    ['A', 'A', '.', '.', 'H', 'I'],
    ['C', 'C', 'C', '.', 'H', '.'],
    ['.', 'G', 'D', 'D', 'E', 'E'],
    ['.', 'G', '.', 'F', 'F', '.']
]

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

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction < 0:  # Move left
            return car_positions[0][1] + direction >= 0 and board[row][car_positions[0][1] + direction] == '.'
        else:  # Move right
            return car_positions[-1][1] + direction < 6 and board[row][car_positions[-1][1] + direction] == '.'
    else:
        col = car_positions[0][1]
        if direction < 0:  # Move up
            return car_positions[0][0] + direction >= 0 and board[car_positions[0][0] + direction][col] == '.'
        else:  # Move down
            return car_positions[-1][0] + direction < 6 and board[car_positions[-1][0] + direction][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [row[:] for row in board]
    car_length = len(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if is_horizontal(positions):
        row = positions[0][0]
        start_col = positions[0][1] + direction
        for i in range(car_length):
            new_board[row][start_col + i] = car
    else:
        col = positions[0][1]
        start_row = positions[0][0] + direction
        for i in range(car_length):
            new_board[start_row + i][col] = car
            
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    start = initial_board
    queue = deque([(start, [])])
    seen = {board_to_string(start)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) has reached the exit
        if current_board[2][4] == 'A':
            return moves
        
        cars = get_car_positions(current_board)
        
        for car, positions in cars.items():
            is_hor = is_horizontal(positions)
            
            # Try moving left/up
            if can_move(current_board, positions, -1):
                new_board = move_car(current_board, car, positions, -1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    new_moves = moves + [f"{car}-1"]
                    queue.append((new_board, new_moves))
            
            # Try moving right/down
            if can_move(current_board, positions, 1):
                new_board = move_car(current_board, car, positions, 1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    new_moves = moves + [f"{car}+1"]
                    queue.append((new_board, new_moves))

solution = solve_puzzle()
print("<<<" + " ".join(solution) + ">>>")