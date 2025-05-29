from collections import deque
import copy

def create_board():
    board = [
        ['B', 'B', 'F', '.', '.', 'H'],
        ['.', 'E', 'F', '.', '.', 'H'],
        ['.', 'E', 'F', 'A', 'A', 'I'],
        ['.', 'C', 'C', 'G', '.', 'I'],
        ['x', '.', '.', 'G', 'D', 'D'],
        ['.', '.', '.', 'x', '.', '.']
    ]
    return board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, car, row, col, is_vertical, length, delta):
    if is_vertical:
        if delta > 0:  # Moving down
            if row + length - 1 + delta >= len(board):
                return False
            return all(board[row + length + i][col] == '.' for i in range(delta))
        else:  # Moving up
            if row + delta < 0:
                return False
            return all(board[row + delta + i][col] == '.' for i in range(-delta))
    else:
        if delta > 0:  # Moving right
            if col + length - 1 + delta >= len(board[0]):
                return False
            return all(board[row][col + length + i] == '.' for i in range(delta))
        else:  # Moving left
            if col + delta < 0:
                return False
            return all(board[row][col + delta + i] == '.' for i in range(-delta))

def move_car(board, car, row, col, is_vertical, length, delta):
    new_board = copy.deepcopy(board)
    # Clear old position
    for i in range(length):
        if is_vertical:
            new_board[row + i][col] = '.'
        else:
            new_board[row][col + i] = '.'
    # Set new position
    for i in range(length):
        if is_vertical:
            new_board[row + delta + i][col] = car
        else:
            new_board[row][col + delta + i] = car
    return new_board

def find_car(board, car):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == car:
                return i, j
    return None

def get_car_info(board, car):
    row, col = find_car(board, car)
    if row + 1 < len(board) and board[row + 1][col] == car:
        length = sum(1 for i in range(row, len(board)) if i < len(board) and board[i][col] == car)
        return row, col, True, length
    else:
        length = sum(1 for j in range(col, len(board[0])) if j < len(board[0]) and board[row][j] == car)
        return row, col, False, length

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    cars = set('ABCDEFGHI')

    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA reaches right edge)
        row, col = find_car(current_board, 'A')
        if col + 2 == len(current_board[0]):
            return moves

        for car in cars:
            if car not in board_to_string(current_board):
                continue
            row, col, is_vertical, length = get_car_info(current_board, car)
            
            # Try moving the car in both directions
            for delta in range(-5, 6):
                if delta == 0:
                    continue
                if is_valid_move(current_board, car, row, col, is_vertical, length, delta):
                    new_board = move_car(current_board, car, row, col, is_vertical, length, delta)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{car}{'+' if delta > 0 else ''}{delta}"]
                        queue.append((new_board, new_moves))

    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")