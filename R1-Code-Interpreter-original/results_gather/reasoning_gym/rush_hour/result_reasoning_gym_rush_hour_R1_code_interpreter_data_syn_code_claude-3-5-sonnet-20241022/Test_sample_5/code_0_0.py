from collections import deque
import copy

def get_board_state():
    return [
        ['.', '.', 'F', 'G', 'B', 'B'],
        ['.', '.', 'F', 'G', 'C', 'C'],
        ['A', 'A', 'F', 'G', 'H', 'I'],
        ['.', '.', '.', '.', 'H', 'I'],
        ['.', '.', '.', '.', '.', 'x'],
        ['.', 'D', 'D', 'E', 'E', '.']
    ]

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'positions': [], 'orientation': None}
                cars[board[i][j]]['positions'].append((i, j))
    
    for car in cars:
        if cars[car]['positions'][0][0] == cars[car]['positions'][1][0]:
            cars[car]['orientation'] = 'horizontal'
        else:
            cars[car]['orientation'] = 'vertical'
    
    return cars

def can_move(board, car_info, car, direction):
    positions = car_info[car]['positions']
    orientation = car_info[car]['orientation']
    
    if orientation == 'horizontal':
        if direction < 0:  # Move left
            leftmost = min(pos[1] for pos in positions)
            return leftmost > 0 and board[positions[0][0]][leftmost - 1] == '.'
        else:  # Move right
            rightmost = max(pos[1] for pos in positions)
            return rightmost < 5 and board[positions[0][0]][rightmost + 1] == '.'
    else:  # vertical
        if direction < 0:  # Move up
            topmost = min(pos[0] for pos in positions)
            return topmost > 0 and board[topmost - 1][positions[0][1]] == '.'
        else:  # Move down
            bottommost = max(pos[0] for pos in positions)
            return bottommost < 5 and board[bottommost + 1][positions[0][1]] == '.'

def move_car(board, car_info, car, direction):
    new_board = [row[:] for row in board]
    positions = car_info[car]['positions']
    orientation = car_info[car]['orientation']
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if orientation == 'horizontal':
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_solved(board):
    # Check if red car (AA) is at the exit
    return board[2][4] == 'A' and board[2][5] == 'A'

def solve_puzzle():
    initial_board = get_board_state()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        car_info = get_car_info(current_board)
        
        for car in car_info:
            for direction in [-1, 1]:
                if can_move(current_board, car_info, car, direction):
                    new_board = move_car(current_board, car_info, car, direction)
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