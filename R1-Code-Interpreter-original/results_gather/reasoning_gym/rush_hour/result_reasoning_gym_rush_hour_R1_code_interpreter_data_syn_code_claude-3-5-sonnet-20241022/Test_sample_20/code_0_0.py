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

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def can_move(board, car_pos, direction):
    if is_horizontal(car_pos):
        if direction > 0:  # move right
            new_pos = car_pos[-1][1] + 1
            return new_pos < 6 and board[car_pos[0][0]][new_pos] == '.'
        else:  # move left
            new_pos = car_pos[0][1] - 1
            return new_pos >= 0 and board[car_pos[0][0]][new_pos] == '.'
    else:
        if direction > 0:  # move down
            new_pos = car_pos[-1][0] + 1
            return new_pos < 6 and board[new_pos][car_pos[0][1]] == '.'
        else:  # move up
            new_pos = car_pos[0][0] - 1
            return new_pos >= 0 and board[new_pos][car_pos[0][1]] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [list(row) for row in board]
    # Clear current position
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new position
    if is_horizontal(car_pos):
        for i, pos in enumerate(car_pos):
            new_board[pos[0]][pos[1] + direction] = car
    else:
        for i, pos in enumerate(car_pos):
            new_board[pos[0] + direction][pos[1]] = car
    
    return [''.join(row) for row in new_board]

def board_to_string(board):
    return ''.join(board)

def solve_puzzle():
    initial_board = [
        'BB.HJK',
        '..xHJK',
        'AAGIJK',
        '.FGI..',
        '.FCC..',
        '.DDEE.'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        cars = get_car_positions(current_board)
        if any(pos[1] == 5 for pos in cars['A']):
            return moves
        
        # Try all possible moves
        for car, positions in cars.items():
            # Try moving right/down
            if can_move(current_board, positions, 1):
                new_board = move_car(current_board, car, positions, 1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, moves + [f"{car}+1"]))
            
            # Try moving left/up
            if can_move(current_board, positions, -1):
                new_board = move_car(current_board, car, positions, -1)
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, moves + [f"{car}-1"]))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')