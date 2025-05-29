from collections import deque
import copy

def get_board_string(board):
    return '\n'.join(''.join(row) for row in board)

def is_solved(board):
    # Check if AA can reach the exit
    for row in board:
        if 'A' in row and row.index('A') == len(row)-2:
            return True
    return False

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
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
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            # Try moving left
            if left >= 0 and board[row][left] == '.':
                moves.append((car, -1))
            
            # Try moving right
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            # Try moving up
            if top >= 0 and board[top][col] == '.':
                moves.append((car, -1))
            
            # Try moving down
            if bottom < len(board) and board[bottom][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, car, direction):
    new_board = [list(row) for row in board]
    positions = cars[car]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return new_board

def solve_puzzle():
    initial_board = [
        list("BBBKCC"),
        list("DDJK.L"),
        list("I.JAAL"),
        list("IEE.xM"),
        list("FF...M"),
        list("GG.xHH")
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(get_board_string(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        cars = get_car_positions(current_board)
        valid_moves = get_valid_moves(current_board, cars)
        
        for car, direction in valid_moves:
            new_board = apply_move(current_board, cars, car, direction)
            board_string = get_board_string(new_board)
            
            if board_string not in visited:
                visited.add(board_string)
                new_moves = moves + [(car, direction)]
                queue.append((new_board, new_moves))
    
    return None

# Solve and format the output
solution = solve_puzzle()
if solution:
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}1")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("No solution found")