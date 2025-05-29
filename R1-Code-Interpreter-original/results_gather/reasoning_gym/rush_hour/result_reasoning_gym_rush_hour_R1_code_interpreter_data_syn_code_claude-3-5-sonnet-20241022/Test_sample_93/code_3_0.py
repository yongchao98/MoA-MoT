from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def find_red_car(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                return i, j
    return None

def get_path_score(board):
    # Find red car position
    red_pos = find_red_car(board)
    if not red_pos:
        return float('inf')
    
    row, col = red_pos
    score = 0
    # Count blocking vehicles to the right of red car
    for j in range(col + 2, len(board[0])):
        if board[row][j] != '.':
            score += 1
    return score

def get_vehicle_length(board, row, col, vehicle):
    # Horizontal check
    if col + 1 < len(board[0]) and board[row][col + 1] == vehicle:
        length = 2
        if col + 2 < len(board[0]) and board[row][col + 2] == vehicle:
            length = 3
        return length, 'H'
    # Vertical check
    elif row + 1 < len(board) and board[row + 1][col] == vehicle:
        length = 2
        if row + 2 < len(board) and board[row + 2][col] == vehicle:
            length = 3
        return length, 'V'
    return 1, None

def make_move(board, row, col, vehicle, length, orientation, direction):
    new_board = [row[:] for row in board]
    
    # Clear current position
    if orientation == 'H':
        for j in range(length):
            new_board[row][col + j] = '.'
        # Place in new position
        if direction == 'L' and col > 0:
            for j in range(length):
                new_board[row][col + j - 1] = vehicle
            return new_board
        elif direction == 'R' and col + length < len(board[0]):
            for j in range(length):
                new_board[row][col + j + 1] = vehicle
            return new_board
    else:  # vertical
        for i in range(length):
            new_board[row + i][col] = '.'
        # Place in new position
        if direction == 'U' and row > 0:
            for i in range(length):
                new_board[row + i - 1][col] = vehicle
            return new_board
        elif direction == 'D' and row + length < len(board):
            for i in range(length):
                new_board[row + i + 1][col] = vehicle
            return new_board
    return None

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    start = (get_path_score(board), board_to_string(board), [])
    visited = {start[1]}
    queue = [start]
    
    while queue:
        _, current_board_str, moves = queue.pop(0)
        current_board = parse_board(current_board_str)
        
        # Check if solved
        red_pos = find_red_car(current_board)
        if red_pos and current_board[red_pos[0]][-1] == 'A':
            return moves
        
        # Try all possible moves for each vehicle
        for i in range(len(current_board)):
            for j in range(len(current_board[0])):
                if current_board[i][j] not in '.x':
                    vehicle = current_board[i][j]
                    length, orientation = get_vehicle_length(current_board, i, j, vehicle)
                    
                    if orientation == 'H':
                        # Try left
                        if j > 0 and current_board[i][j-1] == '.':
                            new_board = make_move(current_board, i, j, vehicle, length, 'H', 'L')
                            if new_board:
                                new_board_str = board_to_string(new_board)
                                if new_board_str not in visited:
                                    visited.add(new_board_str)
                                    score = get_path_score(new_board)
                                    queue.append((score, new_board_str, moves + [f"{vehicle}-1"]))
                        
                        # Try right
                        if j + length < len(current_board[0]) and current_board[i][j+length] == '.':
                            new_board = make_move(current_board, i, j, vehicle, length, 'H', 'R')
                            if new_board:
                                new_board_str = board_to_string(new_board)
                                if new_board_str not in visited:
                                    visited.add(new_board_str)
                                    score = get_path_score(new_board)
                                    queue.append((score, new_board_str, moves + [f"{vehicle}+1"]))
                    
                    elif orientation == 'V':
                        # Try up
                        if i > 0 and current_board[i-1][j] == '.':
                            new_board = make_move(current_board, i, j, vehicle, length, 'V', 'U')
                            if new_board:
                                new_board_str = board_to_string(new_board)
                                if new_board_str not in visited:
                                    visited.add(new_board_str)
                                    score = get_path_score(new_board)
                                    queue.append((score, new_board_str, moves + [f"{vehicle}-1"]))
                        
                        # Try down
                        if i + length < len(current_board) and current_board[i+length][j] == '.':
                            new_board = make_move(current_board, i, j, vehicle, length, 'V', 'D')
                            if new_board:
                                new_board_str = board_to_string(new_board)
                                if new_board_str not in visited:
                                    visited.add(new_board_str)
                                    score = get_path_score(new_board)
                                    queue.append((score, new_board_str, moves + [f"{vehicle}+1"]))
                    
                    # Skip to next vehicle
                    if orientation == 'H':
                        j += length - 1
        
        queue.sort()  # Sort by path score
    
    return None

initial_board = """.BBICC
xDDI.J
x.HAAJ
G.H..K
GEE..K
...FFF"""

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")