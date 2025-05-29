from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        # Determine if car is horizontal or vertical
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Try moving right
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Try moving down
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, car, direction):
    new_board = [list(row) for row in board]
    positions = cars[car]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        '....BB',
        '...x.F',
        'AA.E.F',
        '...E..',
        'CCDD..',
        '......'
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (AA is at the exit)
        if 'A' in cars and any(pos[1] == len(current_board[0])-1 for pos in cars['A']):
            return moves
        
        # Get and apply valid moves
        valid_moves = get_valid_moves(current_board, cars)
        for car, direction in valid_moves:
            new_board = apply_move(current_board, cars, car, direction)
            board_tuple = tuple(new_board)
            
            if board_tuple not in visited:
                visited.add(board_tuple)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format output
solution = solve_puzzle()
if solution:
    formatted_moves = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves += f'{car}{sign}1 '
    print(f'<<<{formatted_moves.strip()}>>>')
else:
    print('<<<No solution found>>>')