from collections import deque
import copy

def create_board():
    return [
        ['B','B','B','C','C','.'],
        ['.','.','H','I','D','D'],
        ['A','A','H','I','J','K'],
        ['G','E','E','E','J','K'],
        ['G','.','F','F','.','.'],
        ['.','.','.','.','.','.']
    ]

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i,j))
    return {k: sorted(v) for k, v in cars.items()}

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, car_pos, direction):
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in car_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in car_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_pos[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in car_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in car_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    # Set new positions
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        cols = [pos[1] for pos in car_pos]
        for col in cols:
            new_board[row][col + direction] = car
    else:
        col = car_pos[0][1]
        rows = [pos[0] for pos in car_pos]
        for row in rows:
            new_board[row + direction][col] = car
    return new_board

def get_blocking_cars(board, red_row, red_right):
    blocking = []
    for col in range(red_right + 1, 6):
        if board[red_row][col] != '.' and board[red_row][col] not in blocking:
            blocking.append(board[red_row][col])
    return blocking

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if red car can reach exit
        if 'A' in cars:
            red_car = cars['A']
            red_row = red_car[0][0]
            red_right = max(pos[1] for pos in red_car)
            
            # If red car is at position 4, we've won
            if red_right == 4:
                return moves
            
            # Check if path to exit is clear
            path_clear = True
            for col in range(red_right + 1, 6):
                if current_board[red_row][col] != '.':
                    path_clear = False
                    break
            
            if path_clear and red_right < 4:
                return moves + [f"A+{4-red_right}"]
        
        # Try all possible moves
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        
                        # Prioritize moves of blocking cars
                        if car in get_blocking_cars(current_board, 2, red_right):
                            queue.appendleft((new_board, new_moves))
                        else:
                            queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")