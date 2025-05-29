from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_valid_move(board, car_pos, direction, distance):
    orientation = 'H' if car_pos[0][0] == car_pos[1][0] else 'V'
    
    if orientation == 'H':
        row = car_pos[0][0]
        if direction > 0:  # moving right
            for pos in car_pos:
                col = pos[1]
                for d in range(1, distance + 1):
                    if col + d >= 6 or board[row][col + d] != '.':
                        return False
        else:  # moving left
            for pos in car_pos:
                col = pos[1]
                for d in range(1, abs(distance) + 1):
                    if col - d < 0 or board[row][col - d] != '.':
                        return False
    else:  # vertical
        col = car_pos[0][1]
        if direction > 0:  # moving down
            for pos in car_pos:
                row = pos[0]
                for d in range(1, distance + 1):
                    if row + d >= 6 or board[row + d][col] != '.':
                        return False
        else:  # moving up
            for pos in car_pos:
                row = pos[0]
                for d in range(1, abs(distance) + 1):
                    if row - d < 0 or board[row - d][col] != '.':
                        return False
    return True

def move_car(board, car_pos, direction, distance):
    new_board = [list(row) for row in board]
    car = board[car_pos[0][0]][car_pos[0][1]]
    
    # Clear old positions
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if car_pos[0][0] == car_pos[1][0]:  # horizontal
        row = car_pos[0][0]
        for pos in car_pos:
            col = pos[1] + distance
            new_board[row][col] = car
    else:  # vertical
        col = car_pos[0][1]
        for pos in car_pos:
            row = pos[0] + distance
            new_board[row][col] = car
    
    return [''.join(row) for row in new_board]

def is_solved(board):
    # Check if red car (AA) is at the exit
    return 'A' in board[2][-2:]

def solve_puzzle(initial_board):
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        cars = get_car_positions(current_board)
        for car, positions in cars.items():
            # Try moving the car in both directions
            for direction in [-1, 1]:
                # Try different distances
                for distance in range(1, 4):
                    if is_valid_move(current_board, positions, direction, distance):
                        new_board = move_car(current_board, positions, direction, distance)
                        board_tuple = tuple(new_board)
                        
                        if board_tuple not in visited:
                            visited.add(board_tuple)
                            move_str = f"{car}{'+' if direction > 0 else '-'}{distance}"
                            queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = [
    "GBBIJK",
    "G..IJK",
    "AAHI..",
    "..HCCC",
    "..xDD.",
    "EEEFF."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")