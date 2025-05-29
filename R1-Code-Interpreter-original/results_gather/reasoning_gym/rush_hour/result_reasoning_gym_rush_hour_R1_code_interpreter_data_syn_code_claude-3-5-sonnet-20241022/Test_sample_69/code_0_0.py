from collections import deque
import copy

def create_board():
    board = [
        ['G', 'B', 'B', 'B', 'I', '.'],
        ['G', 'C', 'C', '.', 'I', '.'],
        ['A', 'A', 'H', '.', 'I', 'J'],
        ['.', 'x', 'H', 'D', 'D', 'J'],
        ['E', 'E', 'E', '.', '.', 'J'],
        ['.', '.', 'F', 'F', 'F', '.']
    ]
    return board

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
        row = car_pos[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in car_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in car_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_pos[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in car_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in car_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        for pos in car_pos:
            new_board[row][pos[1] + direction] = car
    else:
        col = car_pos[0][1]
        for pos in car_pos:
            new_board[pos[0] + direction][col] = car
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car AA reaches right edge)
        if any(pos[1] == 4 for pos in cars['A']):
            return moves
        
        # Try all possible moves for each car
        for car, positions in cars.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")