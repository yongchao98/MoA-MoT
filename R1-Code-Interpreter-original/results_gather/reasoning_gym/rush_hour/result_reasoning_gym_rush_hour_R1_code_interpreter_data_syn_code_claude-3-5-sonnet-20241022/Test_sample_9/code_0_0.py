from collections import deque
import copy

def create_board():
    board = [
        ['.', 'F', '.', 'H', 'B', 'B'],
        ['.', 'F', 'x', 'H', '.', 'J'],
        ['A', 'A', 'G', '.', '.', 'J'],
        ['E', '.', 'G', 'C', 'C', 'J'],
        ['E', '.', 'G', 'I', '.', '.'],
        ['E', 'D', 'D', 'I', '.', '.']
    ]
    return board

def get_cars(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'pos': [(i, j)], 'orientation': None}
                else:
                    cars[board[i][j]]['pos'].append((i, j))
    
    for car in cars:
        if cars[car]['pos'][0][0] == cars[car]['pos'][1][0]:
            cars[car]['orientation'] = 'H'  # horizontal
        else:
            cars[car]['orientation'] = 'V'  # vertical
    return cars

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_possible_moves(board, cars):
    moves = []
    for car in cars:
        if cars[car]['orientation'] == 'H':
            # Try moving left
            r, c = cars[car]['pos'][0]
            if c > 0 and board[r][c-1] == '.':
                moves.append((car, -1))
            # Try moving right
            r, c = cars[car]['pos'][-1]
            if c < 5 and board[r][c+1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            r, c = cars[car]['pos'][0]
            if r > 0 and board[r-1][c] == '.':
                moves.append((car, -1))
            # Try moving down
            r, c = cars[car]['pos'][-1]
            if r < 5 and board[r+1][c] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, car, direction):
    new_board = [row[:] for row in board]
    new_cars = copy.deepcopy(cars)
    
    # Clear current position
    for pos in new_cars[car]['pos']:
        new_board[pos[0]][pos[1]] = '.'
    
    # Update position
    if new_cars[car]['orientation'] == 'H':
        for i in range(len(new_cars[car]['pos'])):
            new_pos = (new_cars[car]['pos'][i][0], new_cars[car]['pos'][i][1] + direction)
            new_cars[car]['pos'][i] = new_pos
            new_board[new_pos[0]][new_pos[1]] = car
    else:
        for i in range(len(new_cars[car]['pos'])):
            new_pos = (new_cars[car]['pos'][i][0] + direction, new_cars[car]['pos'][i][1])
            new_cars[car]['pos'][i] = new_pos
            new_board[new_pos[0]][new_pos[1]] = car
    
    return new_board, new_cars

def solve_puzzle():
    initial_board = create_board()
    initial_cars = get_cars(initial_board)
    
    queue = deque([(initial_board, initial_cars, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if red car (AA) can exit
        if current_cars['A']['pos'][1][1] == 4:
            moves.append(('A', 1))
            return moves
        
        possible_moves = get_possible_moves(current_board, current_cars)
        for car, direction in possible_moves:
            new_board, new_cars = apply_move(current_board, current_cars, car, direction)
            board_string = board_to_string(new_board)
            
            if board_string not in seen:
                seen.add(board_string)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_cars, new_moves))
    
    return None

def format_solution(moves):
    if not moves:
        return "No solution found"
    
    result = []
    for car, direction in moves:
        result.append(f"{car}{'+' if direction > 0 else ''}{direction}")
    return ' '.join(result)

solution = solve_puzzle()
print(format_solution(solution))