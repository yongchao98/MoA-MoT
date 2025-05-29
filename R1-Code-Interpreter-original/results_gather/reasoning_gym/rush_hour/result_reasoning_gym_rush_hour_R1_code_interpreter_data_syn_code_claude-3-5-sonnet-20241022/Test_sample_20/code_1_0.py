from heapq import heappush, heappop
import copy

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

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

def heuristic(board):
    cars = get_car_positions(board)
    # Distance of red car (AA) from exit
    red_car = cars['A']
    exit_pos = (2, 5)  # Exit position
    return manhattan_distance(red_car[1], exit_pos)

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
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    if is_horizontal(car_pos):
        for pos in car_pos:
            new_board[pos[0]][pos[1] + direction] = car
    else:
        for pos in car_pos:
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
    
    # Priority queue with (priority, moves_count, board, moves)
    queue = [(heuristic(initial_board), 0, initial_board, [])]
    seen = {board_to_string(initial_board)}
    
    while queue:
        _, moves_count, current_board, moves = heappop(queue)
        
        cars = get_car_positions(current_board)
        if any(pos[1] == 5 for pos in cars['A']):
            return moves
            
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}1"]
                        priority = moves_count + 1 + heuristic(new_board)
                        heappush(queue, (priority, moves_count + 1, new_board, new_moves))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')