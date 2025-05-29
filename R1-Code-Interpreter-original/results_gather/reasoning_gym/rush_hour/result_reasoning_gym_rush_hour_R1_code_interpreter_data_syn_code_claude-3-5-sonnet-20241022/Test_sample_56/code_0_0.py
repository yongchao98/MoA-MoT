from collections import deque
import copy

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
    return positions[0][0] == positions[-1][0]

def can_move(board, car_positions, direction):
    if is_horizontal(car_positions):
        row = car_positions[0][0]
        if direction > 0:  # moving right
            next_col = car_positions[-1][1] + 1
            return (next_col < len(board[0]) and 
                   board[row][next_col] == '.')
        else:  # moving left
            next_col = car_positions[0][1] - 1
            return (next_col >= 0 and 
                   board[row][next_col] == '.')
    else:  # vertical
        col = car_positions[0][1]
        if direction > 0:  # moving down
            next_row = car_positions[-1][0] + 1
            return (next_row < len(board) and 
                   board[next_row][col] == '.')
        else:  # moving up
            next_row = car_positions[0][0] - 1
            return (next_row >= 0 and 
                   board[next_row][col] == '.')

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if is_horizontal(positions):
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    return [''.join(row) for row in new_board], new_positions

def is_solved(board):
    # Check if red car (AA) reaches right edge
    for row in board:
        if 'A' in row and row.rindex('A') == len(row) - 1:
            return True
    return False

def solve_puzzle():
    initial_board = [
        "..HBBB",
        "..HICC",
        ".AAIJK",
        ".GDDJK",
        ".GEEJL",
        "FF...L"
    ]
    
    visited = set()
    queue = deque([(initial_board, [], {})])
    
    while queue:
        current_board, moves, _ = queue.popleft()
        board_str = '\n'.join(current_board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if is_solved(current_board):
            return moves
        
        cars = get_car_positions(current_board)
        
        for car, positions in cars.items():
            is_hor = is_horizontal(positions)
            
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board, new_positions = move_car(current_board, car, positions, direction)
                    move_str = f"{car}{'+' if direction > 0 else '-'}1"
                    queue.append((new_board, moves + [move_str], {}))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")