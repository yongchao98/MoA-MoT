from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_blocking_score(board):
    # Find cars blocking the path of AA to exit
    score = 0
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                row = i
                col = j
                break
    
    # Check all positions to the right of AA
    for j in range(col + 2, len(board[0])):
        if board[row][j] != '.':
            score += 1
    return score

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    start_state = (board, [])
    visited = {str(board)}
    queue = [(get_blocking_score(board), start_state)]
    
    while queue:
        _, (current_board, moves) = queue.pop(0)
        
        # Check if solved
        for row in current_board:
            if row.count('A') == 2 and row[-1] == 'A':
                return moves
        
        # Try all possible moves
        for i in range(len(current_board)):
            for j in range(len(current_board[0])):
                if current_board[i][j] not in '.x':
                    vehicle = current_board[i][j]
                    # Find vehicle orientation
                    horizontal = (j + 1 < len(current_board[0]) and current_board[i][j + 1] == vehicle)
                    
                    if horizontal:
                        # Try moving left
                        if j > 0 and current_board[i][j-1] == '.':
                            new_board = [row[:] for row in current_board]
                            length = 2 if j+1 >= len(current_board[0]) or current_board[i][j+2] != vehicle else 3
                            if move_horizontal(new_board, i, j, vehicle, length, -1):
                                board_str = str(new_board)
                                if board_str not in visited:
                                    visited.add(board_str)
                                    new_moves = moves + [f"{vehicle}-1"]
                                    score = get_blocking_score(new_board)
                                    queue.append((score, (new_board, new_moves)))
                        
                        # Try moving right
                        if j + 2 < len(current_board[0]) and current_board[i][j+2] == '.':
                            new_board = [row[:] for row in current_board]
                            length = 2 if j+1 >= len(current_board[0]) or current_board[i][j+2] != vehicle else 3
                            if move_horizontal(new_board, i, j, vehicle, length, 1):
                                board_str = str(new_board)
                                if board_str not in visited:
                                    visited.add(board_str)
                                    new_moves = moves + [f"{vehicle}+1"]
                                    score = get_blocking_score(new_board)
                                    queue.append((score, (new_board, new_moves)))
                    
                    else:  # vertical
                        # Try moving up
                        if i > 0 and current_board[i-1][j] == '.':
                            new_board = [row[:] for row in current_board]
                            length = 2 if i+1 >= len(current_board) or current_board[i+2][j] != vehicle else 3
                            if move_vertical(new_board, i, j, vehicle, length, -1):
                                board_str = str(new_board)
                                if board_str not in visited:
                                    visited.add(board_str)
                                    new_moves = moves + [f"{vehicle}-1"]
                                    score = get_blocking_score(new_board)
                                    queue.append((score, (new_board, new_moves)))
                        
                        # Try moving down
                        if i + 2 < len(current_board) and current_board[i+2][j] == '.':
                            new_board = [row[:] for row in current_board]
                            length = 2 if i+1 >= len(current_board) or current_board[i+2][j] != vehicle else 3
                            if move_vertical(new_board, i, j, vehicle, length, 1):
                                board_str = str(new_board)
                                if board_str not in visited:
                                    visited.add(board_str)
                                    new_moves = moves + [f"{vehicle}+1"]
                                    score = get_blocking_score(new_board)
                                    queue.append((score, (new_board, new_moves)))
        
        queue.sort()  # Sort by blocking score
    
    return None

def move_horizontal(board, row, col, vehicle, length, direction):
    if direction == -1 and col == 0:
        return False
    if direction == 1 and col + length >= len(board[0]):
        return False
    
    # Clear current position
    for j in range(length):
        board[row][col + j] = '.'
    
    # Move to new position
    for j in range(length):
        board[row][col + j + direction] = vehicle
    
    return True

def move_vertical(board, row, col, vehicle, length, direction):
    if direction == -1 and row == 0:
        return False
    if direction == 1 and row + length >= len(board):
        return False
    
    # Clear current position
    for i in range(length):
        board[row + i][col] = '.'
    
    # Move to new position
    for i in range(length):
        board[row + i + direction][col] = vehicle
    
    return True

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