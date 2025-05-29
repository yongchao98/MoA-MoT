from collections import deque
import copy

def create_board():
    board = [
        ['.', '.', 'B', 'B', '.', '.'],
        ['x', 'C', 'C', '.', 'x', '.'],
        ['A', 'A', 'G', '.', '.', '.'],
        ['F', '.', 'G', 'D', 'D', '.'],
        ['F', 'E', 'E', 'H', '.', '.'],
        ['F', '.', '.', 'H', '.', '.']
    ]
    return board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, car, row, col, is_vertical, length, delta):
    if is_vertical:
        new_row = row + delta
        if new_row < 0 or new_row + length - 1 >= len(board):
            return False
        if delta > 0:
            return board[new_row + length - 1][col] == '.'
        return board[new_row][col] == '.'
    else:
        new_col = col + delta
        if new_col < 0 or new_col + length - 1 >= len(board[0]):
            return False
        if delta > 0:
            return board[row][new_col + length - 1] == '.'
        return board[row][new_col] == '.'

def move_car(board, car, row, col, is_vertical, length, delta):
    new_board = copy.deepcopy(board)
    if is_vertical:
        for i in range(length):
            new_board[row + i][col] = '.'
        for i in range(length):
            new_board[row + delta + i][col] = car
    else:
        for i in range(length):
            new_board[row][col + i] = '.'
        for i in range(length):
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
    length = 1
    is_vertical = False
    
    if row + 1 < len(board) and board[row + 1][col] == car:
        is_vertical = True
        while row + length < len(board) and board[row + length][col] == car:
            length += 1
    else:
        while col + length < len(board[0]) and board[row][col + length] == car:
            length += 1
            
    return row, col, is_vertical, length

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    cars = set(c for row in initial_board for c in row if c.isalpha())
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) reached the exit
        red_row, red_col = find_car(current_board, 'A')
        if red_col + 2 >= len(current_board[0]) - 1:
            return moves
        
        for car in cars:
            row, col, is_vertical, length = get_car_info(current_board, car)
            
            for delta in [-1, 1]:
                if is_valid_move(current_board, car, row, col, is_vertical, length, delta):
                    new_board = move_car(current_board, car, row, col, is_vertical, length, delta)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        move_str = f"{car}{'+' if delta > 0 else '-'}{abs(delta)}"
                        queue.append((new_board, moves + [move_str]))

    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")