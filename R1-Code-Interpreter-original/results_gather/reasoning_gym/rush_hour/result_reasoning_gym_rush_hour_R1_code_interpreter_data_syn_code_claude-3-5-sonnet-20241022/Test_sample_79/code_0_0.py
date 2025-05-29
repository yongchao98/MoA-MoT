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

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Try moving left
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            
            # Try moving right
            if right < 5 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Try moving up
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            
            # Try moving down
            if bottom < 5 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, car, delta):
    new_board = [row[:] for row in board]
    positions = cars[car]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    new_positions = []
    for pos in positions:
        new_row = pos[0] if horizontal else pos[0] + delta
        new_col = pos[1] + delta if horizontal else pos[1]
        new_board[new_row][new_col] = car
        new_positions.append((new_row, new_col))
    
    return new_board

def solve_puzzle():
    initial_board = [
        list('.BBCC.'),
        list('xxDDK.'),
        list('GAAJK.'),
        list('GEEJ..'),
        list('.HIJ..'),
        list('.HIFFF')
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        cars = get_car_positions(current_board)
        if any(pos[0] == 2 and pos[1] == 5 for pos in cars['A']):
            return moves
        
        # Generate and try all possible moves
        valid_moves = get_valid_moves(current_board, cars)
        for car, delta in valid_moves:
            new_board = apply_move(current_board, cars, car, delta)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [(car, delta)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format output
solution = solve_puzzle()
if solution:
    formatted_moves = []
    for car, delta in solution:
        sign = '+' if delta > 0 else '-'
        formatted_moves.append(f"{car}{sign}{abs(delta)}")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")