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
        if direction > 0:
            col = max(pos[1] for pos in car_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:
            col = min(pos[1] for pos in car_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_pos[0][1]
        if direction > 0:
            row = max(pos[0] for pos in car_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:
            row = min(pos[0] for pos in car_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_pos, direction):
    new_board = copy.deepcopy(board)
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        for pos in car_pos:
            new_board[row][pos[1] + direction] = car
    else:
        col = car_pos[0][1]
        for pos in car_pos:
            new_board[pos[0] + direction][col] = car
    return new_board

def is_path_clear(board, cars):
    if 'A' not in cars:
        return False
    red_car = cars['A']
    row = red_car[0][0]
    right_most = max(pos[1] for pos in red_car)
    for col in range(right_most + 1, 6):
        if board[row][col] != '.':
            return False
    return True

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if path is clear for red car
        if is_path_clear(current_board, cars):
            # Move red car to exit
            red_car = cars['A']
            current_pos = max(pos[1] for pos in red_car)
            spaces_to_move = 4 - current_pos
            if spaces_to_move > 0:
                return moves + [f"A+{spaces_to_move}"]
        
        # Try moving each car
        for car, positions in cars.items():
            # Prioritize cars blocking red car's path
            if car != 'A' and positions[0][0] == 2:  # Same row as red car
                for direction in [-1, 1]:
                    if can_move(current_board, positions, direction):
                        new_board = move_car(current_board, car, positions, direction)
                        board_str = board_to_string(new_board)
                        if board_str not in visited:
                            visited.add(board_str)
                            new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                            queue.appendleft((new_board, new_moves))
        
        # Try other moves
        for car, positions in cars.items():
            if car != 'A' and positions[0][0] != 2:  # Not on red car's row
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