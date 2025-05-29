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
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_blocking_cars(board, cars):
    # Get cars directly blocking the red car's path
    red_row = cars['A'][0][0]
    red_right = max(pos[1] for pos in cars['A'])
    blocking = []
    for j in range(red_right + 1, 6):
        if board[red_row][j] != '.' and board[red_row][j] != 'x':
            if board[red_row][j] not in blocking:
                blocking.append(board[red_row][j])
    return blocking

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

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        blocking = get_blocking_cars(current_board, cars)
        
        # Check if solved
        if any(pos[0] == 2 and pos[1] == 4 for pos in cars['A']):
            return moves
        
        # First priority: Move blocking cars
        moved_blocking = False
        for car in blocking:
            positions = cars[car]
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.appendleft((new_board, new_moves))  # Priority to blocking cars
                        moved_blocking = True
        
        if moved_blocking:
            continue
        
        # Second priority: Move cars that might help clear the path
        for car, positions in cars.items():
            if car not in blocking and car != 'A':
                for direction in [-1, 1]:
                    if can_move(current_board, positions, direction):
                        new_board = move_car(current_board, car, positions, direction)
                        board_str = board_to_string(new_board)
                        if board_str not in visited:
                            visited.add(board_str)
                            new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                            queue.append((new_board, new_moves))
        
        # Last priority: Move red car if possible
        if 'A' in cars:
            positions = cars['A']
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"A{'+' if direction > 0 else '-'}{abs(direction)}"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")