def get_board_string(board):
    return '\n'.join(''.join(row) for row in board)

def can_move(board, vehicle, row, col, direction):
    # For horizontal vehicles
    if vehicle == 'A' or vehicle == 'B' or vehicle == 'C':
        if direction > 0:  # right
            return col + direction < len(board[0]) and board[row][col + direction] == '.'
        else:  # left
            return col + direction >= 0 and board[row][col + direction] == '.'
    # For vertical vehicles
    else:
        if direction > 0:  # down
            return row + direction < len(board) and board[row + direction][col] == '.'
        else:  # up
            return row + direction >= 0 and board[row + direction][col] == '.'

def move_vehicle(board, vehicle, row, col, direction):
    new_board = [list(row) for row in board]
    
    # Clear current position(s)
    if vehicle == 'A':
        new_board[row][col] = '.'
        new_board[row][col+1] = '.'
    elif vehicle == 'B':
        new_board[row][col] = '.'
        new_board[row][col+1] = '.'
    elif vehicle == 'C':
        new_board[row][col] = '.'
        new_board[row][col+1] = '.'
        new_board[row][col+2] = '.'
    else:  # vertical vehicles
        new_board[row][col] = '.'
        new_board[row+1][col] = '.'
        if vehicle == 'H':  # H is 3 units tall
            new_board[row+2][col] = '.'
    
    # Set new position(s)
    if vehicle == 'A':
        new_board[row][col + direction] = 'A'
        new_board[row][col + direction + 1] = 'A'
    elif vehicle == 'B':
        new_board[row][col + direction] = 'B'
        new_board[row][col + direction + 1] = 'B'
    elif vehicle == 'C':
        new_board[row][col + direction] = 'C'
        new_board[row][col + direction + 1] = 'C'
        new_board[row][col + direction + 2] = 'C'
    else:  # vertical vehicles
        new_board[row + direction][col] = vehicle
        new_board[row + direction + 1][col] = vehicle
        if vehicle == 'H':  # H is 3 units tall
            new_board[row + direction + 2][col] = vehicle
    
    return new_board

def solve_puzzle():
    board = [
        list("..Ex.."),
        list("..EG.I"),
        list("AAFGHI"),
        list(".DF.H."),
        list(".DBBH."),
        list("..CCC.")
    ]
    
    moves = []
    
    # 1. Move F down (+1)
    if can_move(board, 'F', 2, 2, 1):
        board = move_vehicle(board, 'F', 2, 2, 1)
        moves.append("F+1")
    
    # 2. Move H up (-2)
    if can_move(board, 'H', 2, 4, -1) and can_move(board, 'H', 1, 4, -1):
        board = move_vehicle(board, 'H', 2, 4, -1)
        board = move_vehicle(board, 'H', 1, 4, -1)
        moves.append("H-2")
    
    # 3. Move B left (-1)
    if can_move(board, 'B', 4, 2, -1):
        board = move_vehicle(board, 'B', 4, 2, -1)
        moves.append("B-1")
    
    # 4. Move C left (-1)
    if can_move(board, 'C', 5, 2, -1):
        board = move_vehicle(board, 'C', 5, 2, -1)
        moves.append("C-1")
    
    # 5. Move A right (+4)
    for _ in range(4):
        if can_move(board, 'A', 2, 0, 1):
            board = move_vehicle(board, 'A', 2, 0, 1)
    moves.append("A+4")
    
    if moves == ["F+1", "H-2", "B-1", "C-1", "A+4"]:
        print("<<<" + " ".join(moves) + ">>>")
    else:
        print("No solution found")

solve_puzzle()