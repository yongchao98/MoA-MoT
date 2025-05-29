from collections import deque
import copy

def create_board():
    board = [
        ['B', 'B', 'G', 'I', '.', 'K'],
        ['E', 'F', 'G', 'I', 'x', 'K'],
        ['E', 'F', 'A', 'A', 'J', 'L'],
        ['C', 'C', 'H', '.', 'J', 'L'],
        ['.', '.', 'H', 'D', 'D', '.'],
        ['.', '.', '.', '.', '.', '.']
    ]
    return board

def find_cars(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha() and board[i][j] != 'x':
                car = board[i][j]
                if car not in cars:
                    cars[car] = []
                cars[car].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = positions[0][1]
            right = positions[-1][1]
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if right < 5 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = positions[0][0]
            bottom = positions[-1][0]
            
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Try moving down
            if bottom < 5 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [row[:] for row in board]
    new_cars = copy.deepcopy(cars)
    
    horizontal = is_horizontal(positions)
    if horizontal:
        row = positions[0][0]
        for col in range(positions[0][1], positions[-1][1] + 1):
            new_board[row][col] = '.'
        new_positions = [(row, col + direction) for col in range(positions[0][1], positions[-1][1] + 1)]
        for r, c in new_positions:
            new_board[r][c] = car
    else:
        col = positions[0][1]
        for row in range(positions[0][0], positions[-1][0] + 1):
            new_board[row][col] = '.'
        new_positions = [(row + direction, col) for row in range(positions[0][0], positions[-1][0] + 1)]
        for r, c in new_positions:
            new_board[r][c] = car
    
    new_cars[car] = new_positions
    return new_board, new_cars

def solve_puzzle():
    initial_board = create_board()
    initial_cars = find_cars(initial_board)
    
    queue = deque([(initial_board, initial_cars, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        
        # Check if red car (AA) can exit
        if 'A' in current_cars and current_cars['A'][1][1] == 4:
            return moves
        
        valid_moves = get_valid_moves(current_board, current_cars)
        for move in valid_moves:
            new_board, new_cars = apply_move(current_board, current_cars, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(move[0], move[1])]
                queue.append((new_board, new_cars, new_moves))
    
    return None

# Solve the puzzle and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for car, spaces in solution:
        sign = '+' if spaces > 0 else ''
        formatted_solution += f'{car}{sign}{spaces} '
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")