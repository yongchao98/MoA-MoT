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
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, car_pos, direction):
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        if direction > 0:  # moving right
            col = max(pos[1] for pos in car_pos) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # moving left
            col = min(pos[1] for pos in car_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_pos[0][1]
        if direction > 0:  # moving down
            row = max(pos[0] for pos in car_pos) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # moving up
            row = min(pos[0] for pos in car_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [list(row) for row in board]
    # Clear current position
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new position
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
    initial_board = [
        list('x.GBBJ'),
        list('CCGHIJ'),
        list('FAAHIK'),
        list('FDD.IK'),
        list('EE.x..'),
        list('......')
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) can exit
        cars = get_car_positions(current_board)
        if 'A' in cars:
            aa_pos = cars['A']
            if aa_pos[0][0] == 2 and max(pos[1] for pos in aa_pos) == len(current_board[0])-1:
                return moves
        
        # Try moving each car
        for car, positions in cars.items():
            # Try moving left/up (-1) and right/down (+1)
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
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